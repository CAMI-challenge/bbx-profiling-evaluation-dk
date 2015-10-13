#This will compute the Unifrac distance from two CAMI formatted taxonomic profiling inputs. The branch lengths are currently all 1 (between each taxonomic rank)
# example call: /local/cluster/bin/python ComputeUnifracDistance.py -A /nfs1/Koslicki_Lab/koslickd/CommonKmers/Validation/Output/GrinderOnDatabase/experiment9_CommonKmers_classification.txt -B /nfs1/Koslicki_Lab/koslickd/CommonKmers/Validation/Output/GrinderOnDatabase/experiment9_grinder_classification.txt
import numpy
import networkx
import scipy.io as sio
import subprocess
import getopt
import sys

def main(argv):
	input_file1 = ''
	input_file2 = ''
	try:
		opts, args = getopt.getopt(argv,"g:r:h",["GroundTruth=","Reconstruction="])
	except getopt.GetoptError:
		print 'Files must be in the Bioboxes profiling format: https://github.com/bioboxes/rfc/tree/master/data-format \n Call using: python EMDUnifrac.py -g <GroundTruth.cami> -r <Reconstruction.cami>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Files must be in the Bioboxes profiling format: https://github.com/bioboxes/rfc/tree/master/data-format \n Call using: python EMDUnifrac.py -g <GroundTruth.cami> -r <Reconstruction.cami>'
			sys.exit(2)
		elif opt in ("-g", "--GroundTruth"):
			input_file1 = arg
		elif opt in ("-r", "--Reconstruction"):
			input_file2 = arg

	#Read in classification 1
	tax_path1 = list();
	tax_ids1 = list();
	weights1 = dict();
	fid = open(input_file1,'r');
	temp = fid.readlines();
	for line in temp:
		temp_split = line.split();
		if len(temp_split)>1: #skip blank lines
			if '@' not in temp_split[0] and '#' not in temp_split[0]: #not comment or header
				tax_path1.append(temp_split[2]) #add the whole taxpath
				tax_ids1.append(temp_split[0]) #just terminal tax ID
				weights1[temp_split[0]] = float(temp_split[-1]) #the associated weight
	fid.close()

	#Read in classification 2
	tax_path2 = list();
	tax_ids2 = list();
	weights2 = dict();
	fid = open(input_file2,'r');
	temp = fid.readlines();
	for line in temp:
		temp_split = line.split();
		if len(temp_split)>1:
			if '@' not in temp_split[0] and '#' not in temp_split[0]: #not comment or header
				tax_path2.append(temp_split[2]) #add the whole taxpath
				tax_ids2.append(temp_split[0]) #just terminal tax ID
				weights2[temp_split[0]] = float(temp_split[-1]) #the associated weight
	fid.close()
	all_taxpaths = list(set(tax_path1) | set(tax_path2));
	
	# form the graph
	G = networkx.Graph() #undirected for the proper distance matrix
	G.add_node('-1') #root
	Gdir = networkx.DiGraph() #directed so we can keep track of the descendants
	Gdir.add_node('-1') #root
	for tax_path in all_taxpaths:
		tax_ids = tax_path.split('|')
		for i in range(len(tax_ids)-1,-1,-1): #iterate backwards through list
			if tax_ids[i] != '':
				branch_len = 1; #This is where we could add the actual branch lengths
				if i==0: #if it's the first entry, connect it to the root
					parent_id = '-1'
					G.add_edge(parent_id,tax_ids[i], weight = branch_len)
					Gdir.add_edge(parent_id,tax_ids[i], weight = branch_len)
				else: #if it's not the first entry, look for its parent
					for j in range(i-1,-1,-1): #traverse up the taxpath until a non '||' is found
						if tax_ids[j] != '': #found the parent, so add the clade
							parent_id = tax_ids[j]
							G.add_edge(parent_id, tax_ids[i], weight = branch_len)
							Gdir.add_edge(parent_id, tax_ids[i], weight = branch_len)
							break #you found the parent, so stop traversing up the taxpath
						else:
							branch_len = branch_len + 1 #you found a '||' so add one to the taxonomy
	nodes = G.nodes()
	D = numpy.zeros([len(nodes),len(nodes)],dtype = numpy.int8)
	
	#Compute all pairwise distances
	dist_dict = networkx.all_pairs_dijkstra_path_length(G)
	lookup = {};
	for i, elem in enumerate(nodes):
		lookup[elem] = i
	#Turn the dict of dicts into a distance matrix
	for node1 in nodes:
		for node2 in nodes:
			if node2 in dist_dict[node1]:
				D[lookup[node1],lookup[node2]] = dist_dict[node1][node2]
			elif node1 in dist_dict[node]:
				D[lookup[node1],lookup[node2]] = dist_dict[node2][node1]
	#test if symmetric
	#(D.transpose() == D).all()
	
	#Form the probability distribution for input1. The CAMI format automatically sums up over lower taxa, so I need to subtract descendant weights
	prob_dist1 = numpy.zeros(len(nodes), dtype = numpy.float64)
	for i in range(0,len(nodes)):
		if nodes[i] in weights1:
			prob_dist1[i] = weights1[nodes[i]]
			#now subtract all of the descendant weights. This is so we end up with a probability distribution.
			for desc_node in Gdir.neighbors(nodes[i]): #just the one-step descendants
				if desc_node in weights1:
					prob_dist1[i] = max([prob_dist1[i]-weights1[desc_node], 0]) #threshold to zero to get around things like -1.0e-16
	prob_dist1 = prob_dist1/numpy.sum(prob_dist1) #normalize to make sure it's a prob. distribution
	
	#Same thing for second file
	prob_dist2 = numpy.zeros(len(nodes), dtype = numpy.float64)
	for i in range(0,len(nodes)):
		if nodes[i] in weights2:
			prob_dist2[i] = weights2[nodes[i]]
			#now subtract all of the descendant weights. This is so we end up with a probability distribution.
			for desc_node in Gdir.neighbors(nodes[i]): #just the one-step descendants
				if desc_node in weights2:
					prob_dist2[i] = max([prob_dist2[i]-weights2[desc_node], 0]) #threshold to zero to get around things like -1.0e-16
	prob_dist2 = prob_dist2/numpy.sum(prob_dist2) #normalize to make sure it's a prob. distribution
	
	#Now compute the EMD
	
	Distance = D
	MassA = prob_dist1; # Pick off the masses for the support of A and B (We'll be editting these, so we leave A, B alone)
	MassB = prob_dist2;
	EMD = 0; # Initialize the EMD computation.
	d = 0; # Initialize distance to move mass at 0.
	InitialMove = [min([MassA[i], MassB[i]]) for i in range(0,len(MassA))]
	Flow = numpy.zeros([len(MassA),len(MassA)],dtype = numpy.float64) #Later, can make this a sparse matrix for computational efficiency
	MassA = MassA - InitialMove; # Removes mass that was moved initially.
	MassB = MassB - InitialMove;
	
	while sum(MassA) > 1e-10: # While we still have mass to move,
		d = d+1; # increment distance to move
		IndicesSortedMassSourceA = numpy.flipud(numpy.argsort(MassA)); # sort the sources, big to small.
		for SourceA in IndicesSortedMassSourceA: # Now, for each source of mass in A
			if MassA[SourceA]==0: # Have we gotten through all the sources with mass left? If so, break.
				break
			dNeighborsSourceA = numpy.argwhere(Distance[SourceA,:] == d);# Find the n-mers in B which are distance d from our source.
			dNeighborsSourceA = dNeighborsSourceA.flatten()
			IndicesSortedMassB = numpy.flipud(numpy.argsort(MassB[dNeighborsSourceA])); # We order the sinks, we're going to fill from the top down.
			for SinkB in IndicesSortedMassB: # Iterating over the indices of the sinks.
				if MassB[dNeighborsSourceA[SinkB]] == 0:
					break
				if MassA[SourceA] - MassB[dNeighborsSourceA[SinkB]] >= 0: # check to see if our source can fulfil this sink. 
					Flow[SourceA,dNeighborsSourceA[SinkB]] = MassB[dNeighborsSourceA[SinkB]]; # If so, note the flow.
					EMD = EMD + d*MassB[dNeighborsSourceA[SinkB]]; # update the EMD calc,
					MassA[SourceA] = MassA[SourceA] - MassB[dNeighborsSourceA[SinkB]]; # then remove mass from A
					MassB[dNeighborsSourceA[SinkB]] = 0; # and remove mass from B.
				else: # otherwise, we've run out of mass from SourceA.
					Flow[SourceA,dNeighborsSourceA[SinkB]] = MassA[SourceA]; # If so, note the flow to this last sink,
					EMD = EMD + d*MassA[SourceA]; # update the EMD calc,
					MassB[dNeighborsSourceA[SinkB]]=MassB[dNeighborsSourceA[SinkB]]-MassA[SourceA]; # remove mass from B,
					MassA[SourceA] = 0; # then remove mass from A
					break # and end the loop, with no more mass to distribute from this source.
	print EMD

if __name__ == "__main__":
	main(sys.argv[1:])

