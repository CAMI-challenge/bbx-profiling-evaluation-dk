#!/usr/bin/env python

"""
# This will compute the Unifrac distance from two CAMI formatted taxonomic profiling inputs.
The branch lengths are currently all 1 (between each taxonomic rank)
# example call:
/local/cluster/bin/python ComputeUnifracDistance.py
-A experiment9_CommonKmers_classification.txt
-B experiment9_grinder_classification.txt
"""

import numpy
import networkx
import getopt
import sys


def main(argv):
	input_file1 = ''
	input_file2 = ''
	try:
		opts, args = getopt.getopt(argv, "g:r:h", ["GroundTruth=", "Reconstruction="])
	except getopt.GetoptError:
		print 'Files must be in the Bioboxes profiling format: https://github.com/bioboxes/rfc/tree/master/data-format'
		print 'Call using: python EMDUnifrac.py -g <GroundTruth.cami> -r <Reconstruction.cami>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Files must be in the Bioboxes profiling format: https://github.com/bioboxes/rfc/tree/master/data-format'
			print 'Call using: python EMDUnifrac.py -g <GroundTruth.cami> -r <Reconstruction.cami>'
			sys.exit(2)
		elif opt in ("-g", "--GroundTruth"):
			input_file1 = arg
		elif opt in ("-r", "--Reconstruction"):
			input_file2 = arg
	print emd_unifrac(input_file1, input_file2)


def read_taxonomy_file(file_path):
	tax_path = list()
	tax_ids = list()
	weights = dict()
	with open(file_path, 'r') as read_handler:
		for line in read_handler:
			line = line.rstrip()
			if len(line) == 0:
				continue  # skip blank lines
			if line[0] in ['@', '#']:
				continue  # skip comment or header
			temp_split = line.split('\t')
			tax_path.append(temp_split[2])  # add the whole taxpath
			tax_ids.append(temp_split[0])  # just terminal tax ID
			weights[temp_split[0]] = float(temp_split[4])  # the associated weight
	return tax_ids, tax_path, weights


def get_probability_distribution(nodes, network_graph_directed, weights):
	# Form the probability distribution for input1.
	# The CAMI format automatically sums up over lower taxa, so I need to subtract descendant weights
	probability_distribution = numpy.zeros(len(nodes), dtype=numpy.float64)
	for i in range(0, len(nodes)):
		if nodes[i] in weights:
			probability_distribution[i] = weights[nodes[i]]
			# now subtract all of the descendant weights. This is so we end up with a probability distribution.
			for desc_node in network_graph_directed.neighbors(nodes[i]):  # just the one-step descendants
				if desc_node in weights:
					# threshold to zero to get around things like -1.0e-16
					probability_distribution[i] = max([probability_distribution[i]-weights[desc_node], 0])
	# normalize to make sure it's a prob. distribution
	return probability_distribution/numpy.sum(probability_distribution)


def emd_unifrac(file_path_truth, file_path_query):
	# Read in classification 1
	tax_ids_truth, tax_path_truth, weights_truth = read_taxonomy_file(file_path_truth)
	tax_ids_query, tax_path_query, weights_query = read_taxonomy_file(file_path_query)

	all_taxpaths = list(set(tax_path_truth) | set(tax_path_query))
	
	# form the graph
	network_graph = networkx.Graph()  # undirected for the proper distance matrix
	network_graph_directed = networkx.DiGraph()  # directed so we can keep track of the descendants
	for tax_path in all_taxpaths:
		tax_ids = tax_path.split('|')
		for i in range(len(tax_ids) - 1, -1, -1):  # iterate backwards through list
			if tax_ids[i] == '':
				continue
			parent_id = '-1'  # in case of unknown path
			branch_len = 1  # This is where we could add the actual branch lengths
			for j in range(i-1, -1, -1):  # traverse up the taxpath until a non '||' is found
				if tax_ids[j] == '':
					branch_len += 1  # you found a '||' so add one to the taxonomy
					continue
				# found the parent, so add the clade
				parent_id = tax_ids[j]
				break
			network_graph.add_edge(parent_id, tax_ids[i], weight=branch_len)
			network_graph_directed.add_edge(parent_id, tax_ids[i], weight=branch_len)
			# break  # you found the parent, so stop traversing up the taxpath
	nodes = network_graph.nodes()
	node_distance_matrix = numpy.zeros([len(nodes), len(nodes)], dtype=numpy.int8)

	# Compute all pairwise distances
	dist_dict = networkx.all_pairs_dijkstra_path_length(network_graph)
	lookup = {}
	for i, elem in enumerate(nodes):
		lookup[elem] = i
	# Turn the dict of dicts into a distance matrix
	for node1 in nodes:
		for node2 in nodes:
			if node2 in dist_dict[node1]:
				node_distance_matrix[lookup[node1], lookup[node2]] = dist_dict[node1][node2]
			elif node1 in dist_dict[node2]:
				node_distance_matrix[lookup[node1], lookup[node2]] = dist_dict[node2][node1]
	# test if symmetric
	# (D.transpose() == D).all()
	
	# Form the probability distribution for input1.
	# The CAMI format automatically sums up over lower taxa, so I need to subtract descendant weights
	prob_dist_truth = get_probability_distribution(nodes, network_graph_directed, weights_truth)
	#  Same thing for second file
	prob_dist_query = get_probability_distribution(nodes, network_graph_directed, weights_query)
	return compute_emd(prob_dist_truth, prob_dist_query, node_distance_matrix)


def compute_emd(mass_truth, mass_query, node_distance_matrix):
	# Now compute the EMD
	# Pick off the masses for the support of A and B (We'll be editting these, so we leave A, B alone) # why leave alone??
	# removing mass at same location in both (distance 0)
	initial_move = [min([mass_truth[i], mass_query[i]]) for i in range(0, len(mass_truth))]
	# Later, can make this a sparse matrix for computational efficiency
	mass_truth = mass_truth - initial_move  # Removes mass that was moved initially.
	mass_query = mass_query - initial_move

	flow = numpy.zeros([len(mass_truth), len(mass_truth)], dtype=numpy.float64)
	emd = 0  # Initialize the EMD computation.
	d = 0  # Initialize distance to move mass at 0.
	max_distance = numpy.max(node_distance_matrix)
	# while sum(mass_truth) > 1e-10:  # While we still have mass to move,
	while d < max_distance:  # While we still have mass to move,
		d += 1  # increment distance to move
		indices_sorted_mass_source_truth = numpy.flipud(numpy.argsort(mass_truth))  # sort the sources, big to small.
		for indice_source_truth in indices_sorted_mass_source_truth:  # Now, for each source of mass in A
			if mass_truth[indice_source_truth] == 0:  # Have we gotten through all the sources with mass left? If so, break.
				break
			# Find the n-mers in B which are distance d from our source.
			# d_neighbors_source_a = numpy.argwhere(Distance[source_a, :] == d)
			d_neighbors_truth = numpy.argwhere(node_distance_matrix[indice_source_truth, :] == d)
			d_neighbors_truth = d_neighbors_truth.flatten()
			# We order the sinks, we're going to fill from the top down.
			indices_sorted_mass_b = numpy.flipud(numpy.argsort(mass_query[d_neighbors_truth]))
			for sink_b in indices_sorted_mass_b:  # Iterating over the indices of the sinks.
				if mass_query[d_neighbors_truth[sink_b]] == 0:
					break
				# check to see if our source can fulfil this sink.
				if mass_truth[indice_source_truth] - mass_query[d_neighbors_truth[sink_b]] >= 0:
					# If so, note the flow.
					flow[indice_source_truth, d_neighbors_truth[sink_b]] = mass_query[d_neighbors_truth[sink_b]]
					emd += d*mass_query[d_neighbors_truth[sink_b]]  # update the EMD calc,
					# then remove mass from A
					mass_truth[indice_source_truth] = mass_truth[indice_source_truth] - mass_query[d_neighbors_truth[sink_b]]
					mass_query[d_neighbors_truth[sink_b]] = 0  # and remove mass from B.
				else:  # otherwise, we've run out of mass from SourceA.
					# If so, note the flow to this last sink,
					flow[indice_source_truth, d_neighbors_truth[sink_b]] = mass_truth[indice_source_truth]
					emd += d*mass_truth[indice_source_truth]  # update the EMD calc,
					# remove mass from B,
					mass_query[d_neighbors_truth[sink_b]] = mass_query[d_neighbors_truth[sink_b]] - mass_truth[indice_source_truth]
					mass_truth[indice_source_truth] = 0  # then remove mass from A
					break  # and end the loop, with no more mass to distribute from this source.
	if numpy.sum(mass_truth) > 1e-10:
		raise ArithmeticError("Mass moving failed: {}".format(numpy.sum(mass_truth)))
	return emd


if __name__ == "__main__":
	main(sys.argv[1:])
