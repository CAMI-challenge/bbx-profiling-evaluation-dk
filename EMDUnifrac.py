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


def emd_unifrac(input_file_truth, input_file_query):
	# Read in classification 1
	tax_path_truth = list()
	tax_ids_truth = list()
	weights_truth = dict()
	with open(input_file_truth, 'r') as fid:
		for line in fid:
			line = line.rstrip()
			if len(line) == 0:
				continue  # skip blank lines
			if line[0] in ['@', '#']:
				continue  # skip comment or header
			temp_split = line.split('\t')
			tax_path_truth.append(temp_split[2])  # add the whole taxpath
			tax_ids_truth.append(temp_split[0])  # just terminal tax ID
			weights_truth[temp_split[0]] = float(temp_split[4])  # the associated weight

	# Read in classification 2
	tax_path_query = list()
	tax_ids_query = list()
	weights_query = dict()
	with open(input_file_query, 'r') as fid:
		for line in fid:
			line = line.rstrip()
			if len(line) == 0:
				continue  # skip blank lines
			if line[0] in ['@', '#']:
				continue  # skip comment or header
			temp_split = line.split('\t')
			tax_path_query.append(temp_split[2])  # add the whole taxpath
			tax_ids_query.append(temp_split[0])  # just terminal tax ID
			weights_query[temp_split[0]] = float(temp_split[4])  # the associated weight

	all_taxpaths = list(set(tax_path_truth) | set(tax_path_query))
	
	# form the graph
	network_graph = networkx.Graph()  # undirected for the proper distance matrix
	network_graph.add_node('-1')  # root
	network_graph_directed = networkx.DiGraph()  # directed so we can keep track of the descendants
	network_graph_directed.add_node('-1')  # root
	for tax_path in all_taxpaths:
		tax_ids = tax_path.split('|')
		for i in range(len(tax_ids) - 1, -1, -1):  # iterate backwards through list
			if tax_ids[i] == '':
				continue
			branch_len = 1  # This is where we could add the actual branch lengths
			if i == 0:  # if it's the first entry, connect it to the root
				parent_id = '-1'
				network_graph.add_edge(parent_id, tax_ids[i], weight=branch_len)
				network_graph_directed.add_edge(parent_id, tax_ids[i], weight=branch_len)
			else:  # if it's not the first entry, look for its parent
				for j in range(i-1, -1, -1):  # traverse up the taxpath until a non '||' is found
					if tax_ids[j] == '':  # found the parent, so add the clade
						branch_len += 1  # you found a '||' so add one to the taxonomy
						continue
					parent_id = tax_ids[j]
					network_graph.add_edge(parent_id, tax_ids[i], weight=branch_len)
					network_graph_directed.add_edge(parent_id, tax_ids[i], weight=branch_len)
					break  # you found the parent, so stop traversing up the taxpath
	nodes = network_graph.nodes()
	tax_distance_matrix = numpy.zeros([len(nodes), len(nodes)], dtype=numpy.int8)

	# Compute all pairwise distances
	dist_dict = networkx.all_pairs_dijkstra_path_length(network_graph)
	lookup = {}
	for i, elem in enumerate(nodes):
		lookup[elem] = i
	# Turn the dict of dicts into a distance matrix
	for node1 in nodes:
		for node2 in nodes:
			if node2 in dist_dict[node1]:
				tax_distance_matrix[lookup[node1], lookup[node2]] = dist_dict[node1][node2]
			elif node1 in dist_dict[node2]:
				tax_distance_matrix[lookup[node1], lookup[node2]] = dist_dict[node2][node1]
	# test if symmetric
	# (D.transpose() == D).all()
	
	# Form the probability distribution for input1.
	# The CAMI format automatically sums up over lower taxa, so I need to subtract descendant weights
	prob_dist_truth = numpy.zeros(len(nodes), dtype=numpy.float64)
	for i in range(0, len(nodes)):
		if nodes[i] in weights_truth:
			prob_dist_truth[i] = weights_truth[nodes[i]]
			# now subtract all of the descendant weights. This is so we end up with a probability distribution.
			for desc_node in network_graph_directed.neighbors(nodes[i]):  # just the one-step descendants
				if desc_node in weights_truth:
					# threshold to zero to get around things like -1.0e-16
					prob_dist_truth[i] = max([prob_dist_truth[i]-weights_truth[desc_node], 0])
	prob_dist_truth = prob_dist_truth/numpy.sum(prob_dist_truth)  # normalize to make sure it's a prob. distribution
	
	#  Same thing for second file
	prob_dist_query = numpy.zeros(len(nodes), dtype=numpy.float64)
	for i in range(0, len(nodes)):
		if nodes[i] in weights_query:
			prob_dist_query[i] = weights_query[nodes[i]]
			#  now subtract all of the descendant weights. This is so we end up with a probability distribution.
			for desc_node in network_graph_directed.neighbors(nodes[i]):  # just the one-step descendants
				if desc_node in weights_query:
					# threshold to zero to get around things like -1.0e-16
					prob_dist_query[i] = max([prob_dist_query[i]-weights_query[desc_node], 0])
	prob_dist_query = prob_dist_query/numpy.sum(prob_dist_query)  # normalize to make sure it's a prob. distribution
	
	# Now compute the EMD
	# Pick off the masses for the support of A and B (We'll be editting these, so we leave A, B alone)
	# Distance = distance_matrix
	mass_truth = prob_dist_truth
	mass_query = prob_dist_query
	initial_move = [min([mass_truth[i], mass_query[i]]) for i in range(0, len(mass_truth))]
	# Later, can make this a sparse matrix for computational efficiency
	flow = numpy.zeros([len(mass_truth), len(mass_truth)], dtype=numpy.float64)
	mass_truth = mass_truth - initial_move  # Removes mass that was moved initially.
	mass_query = mass_query - initial_move

	emd = 0  # Initialize the EMD computation.
	d = 0  # Initialize distance to move mass at 0.
	while sum(mass_truth) > 1e-10:  # While we still have mass to move,
		d += 1  # increment distance to move
		indices_sorted_mass_source_truth = numpy.flipud(numpy.argsort(mass_truth))  # sort the sources, big to small.
		for indice_source_truth in indices_sorted_mass_source_truth:  # Now, for each source of mass in A
			if mass_truth[indice_source_truth] == 0:  # Have we gotten through all the sources with mass left? If so, break.
				break
			# Find the n-mers in B which are distance d from our source.
			# d_neighbors_source_a = numpy.argwhere(Distance[source_a, :] == d)
			d_neighbors_source_a = numpy.argwhere(tax_distance_matrix[indice_source_truth, :] == d)
			d_neighbors_source_a = d_neighbors_source_a.flatten()
			# We order the sinks, we're going to fill from the top down.
			indices_sorted_mass_b = numpy.flipud(numpy.argsort(mass_query[d_neighbors_source_a]))
			for sink_b in indices_sorted_mass_b:  # Iterating over the indices of the sinks.
				if mass_query[d_neighbors_source_a[sink_b]] == 0:
					break
				# check to see if our source can fulfil this sink.
				if mass_truth[indice_source_truth] - mass_query[d_neighbors_source_a[sink_b]] >= 0:
					# If so, note the flow.
					flow[indice_source_truth, d_neighbors_source_a[sink_b]] = mass_query[d_neighbors_source_a[sink_b]]
					emd += d*mass_query[d_neighbors_source_a[sink_b]]  # update the EMD calc,
					# then remove mass from A
					mass_truth[indice_source_truth] = mass_truth[indice_source_truth] - mass_query[d_neighbors_source_a[sink_b]]
					mass_query[d_neighbors_source_a[sink_b]] = 0  # and remove mass from B.
				else:  # otherwise, we've run out of mass from SourceA.
					# If so, note the flow to this last sink,
					flow[indice_source_truth, d_neighbors_source_a[sink_b]] = mass_truth[indice_source_truth]
					emd += d*mass_truth[indice_source_truth]  # update the EMD calc,
					# remove mass from B,
					mass_query[d_neighbors_source_a[sink_b]] = mass_query[d_neighbors_source_a[sink_b]]-mass_truth[indice_source_truth]
					mass_truth[indice_source_truth] = 0  # then remove mass from A
					break  # and end the loop, with no more mass to distribute from this source.
	return emd


if __name__ == "__main__":
	main(sys.argv[1:])
