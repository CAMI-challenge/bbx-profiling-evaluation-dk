import os
import subprocess
import sys
import getopt

#input_file1='../EMD/Test/test_truth_CAMI.txt'
#input_file2='../EMD/Test/test_reconstruction_CAMI.txt'
#output_file='test.txt'
#EMD_loc = "../EMD/EMDUnifrac.py"
#python_loc = "/local/cluster/bin/python"
#epsilon = 0
#normalize="n"

def main(argv):
	input_file1 = ''
	input_file2 = ''
	try:
		opts, args = getopt.getopt(argv,"g:r:o:u:p:e:n:h",["GroundTruth=","Reconstruction=","Output=","UnifracLocation=","PythonLocation=","Epsilon=","Normalize="])
	except getopt.GetoptError:
		print 'Call using: python ProfilingMetrics.py -g <GroundTruth.profile> -r <Reconstruction.profile> -o <Output.txt> -u /path/to/UnifracEMD/script -p /path/to/python -e epsilon -n <normalize(y/n)'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Call using: python ProfilingMetrics.py -g <GroundTruth.profile> -r <Reconstruction.profile> -o <Output.txt> -u /path/to/UnifracEMD/script -p /path/to/python -e epsilon -n <normalize(y/n)'
			sys.exit(2)
		elif opt in ("-g", "--GroundTruth"):
			input_file1 = arg
		elif opt in ("-r", "--Reconstruction"):
			input_file2 = arg
		elif opt in ("-o", "--Output"):
			output_file = arg
		elif opt in ("-u", "--UnifracLocation"):
			EMD_loc = arg
		elif opt in ("-p", "--PythonLocation"):
			python_loc = arg
		elif opt in ("-e", "--Epsilon"):
			epsilon = float(arg)
		elif opt in ("-n", "--Normalize"):
			normalize = arg

	tax_path1 = list();
	tax_ids1 = list();
	weights1_dict = dict();
	fid = open(input_file1,'r');
	temp = fid.readlines();
	for line in temp:
		temp_split = line.split();
		if len(temp_split)>1: #skip blank lines
			if line.lower().split(':')[0]=='@ranks':
				ranks1 = line.strip().split(':')[1].split('|')
			if '@' not in temp_split[0] and '#' not in temp_split[0]: #not comment or header
				tax_path1.append(temp_split[2]) #add the whole taxpath
				tax_ids1.append(temp_split[0]) #just terminal tax ID
				weights1_dict[temp_split[0]] = float(temp_split[-1]) #the associated weight
	
	fid.close()
	#Read in classification 2
	tax_path2 = list();
	tax_ids2 = list();
	weights2_dict = dict();
	fid = open(input_file2,'r');
	temp = fid.readlines();
	for line in temp:
		temp_split = line.split();
		if len(temp_split)>1:
			if line.lower().split(':')[0]=='@ranks':
				ranks2 = line.strip().split(':')[1].split('|')
			if '@' not in temp_split[0] and '#' not in temp_split[0]: #not comment or header
				tax_path2.append(temp_split[2]) #add the whole taxpath
				tax_ids2.append(temp_split[0]) #just terminal tax ID
				weights2_dict[temp_split[0]] = float(temp_split[-1]) #the associated weight
	
	fid.close()
	all_taxpaths = list(set(tax_path1) | set(tax_path2));
	
	#Use cutoff 
	for tax_id in tax_ids1:
		if weights1_dict[tax_id] < epsilon:
			weights1_dict[tax_id] = 0
	
	#Delete zero weighted taxIDs
	reduced_tax_ids1 = list()
	reduced_tax_path1 = list()
	reduced_weights1_dict = dict();
	for key in weights1_dict.keys():
		if weights1_dict[key] !=0:
			reduced_tax_ids1.append(key)
			reduced_tax_path1.append(tax_path1[tax_ids1.index(key)])
			reduced_weights1_dict[key] = weights1_dict[key]
	
	
	#Use cutoff 
	for tax_id in tax_ids2:
		if weights2_dict[tax_id] < epsilon:
			weights2_dict[tax_id] = 0
	
	#Delete zero weighted taxIDs
	reduced_tax_ids2 = list()
	reduced_tax_path2 = list()
	reduced_weights2_dict = dict();
	for key in weights2_dict.keys():
		if weights2_dict[key] !=0:
			reduced_tax_ids2.append(key)
			reduced_tax_path2.append(tax_path2[tax_ids2.index(key)])
			reduced_weights2_dict[key] = weights2_dict[key]
	
	#Get the ranks, whichever is longest
	if len(ranks1) > len(ranks2):
		rank_names = ranks1
	else:
		rank_names = ranks2
	
	
	#open the output file
	fid = open(output_file,'w')
	fid.write("@input\t " + os.path.basename(input_file2) + "\n")
	fid.write("@ground_truth\t " + os.path.basename(input_file1) + "\n")
	fid.write("@Taxanomic ranks\t "+ '|'.join(rank_names) + "\n")
	
	TPs=list()
	FPs=list()
	FNs=list()
	Precs=list()
	Senss=list()
	L1norms=list()
	
	for rank in rank_names:
		#select tax IDs and freqs of interest from input1
		rank_pos = rank_names.index(rank)
		#get the tax ids at the rank of interest
		rank_taxIDs1 = list();
		rank_freqs1 = list();
		for i in range(0,len(reduced_tax_path1)):
			temp = reduced_tax_path1[i].split('|')
			if len(temp)==rank_pos + 1:
				rank_taxIDs1.append(temp[rank_pos])
				rank_freqs1.append(reduced_weights1_dict[temp[rank_pos]])
		#select tax IDs and freqs of interest from input2
		rank_taxIDs2 = list();
		rank_freqs2 = list();
		for i in range(0,len(reduced_tax_path2)):
			temp = reduced_tax_path2[i].split('|')
			if len(temp)==rank_pos + 1:
				rank_taxIDs2.append(temp[rank_pos])
				rank_freqs2.append(reduced_weights2_dict[temp[rank_pos]])
		#normalize if this was asked for
		if normalize == "y":
			if sum(rank_freqs1)>0:
				rank_freqs1 = [rank_freqs1[i] / sum(rank_freqs1) for i in range(0,len(rank_freqs1))]
			if sum(rank_freqs2)>0:
				rank_freqs2 = [rank_freqs2[i] / sum(rank_freqs2) for i in range(0,len(rank_freqs2))]
		set_rank_taxIDs1 = set(rank_taxIDs1)
		set_rank_taxIDs2 = set(rank_taxIDs2)
		if sum(rank_freqs1)+sum(rank_freqs2) > 0:
			#TP/FP/etc
			TP = len(set_rank_taxIDs1 & set_rank_taxIDs2) #in both
			FP = len(set_rank_taxIDs2.difference(set_rank_taxIDs1)) #in reconstruction not truth
			FN = len(set_rank_taxIDs1.difference(set_rank_taxIDs2)) #in truth not reconstruction
			if TP+FP==0:
				Prec = -1
			else:
				Prec = float(TP)/(TP+FP)
			if TP+FN==0:
				Sens = -1
			else:
				Sens = float(TP)/(TP+FN)
			#L1 error
			all_taxIDs = set_rank_taxIDs1 | set_rank_taxIDs2
			L1norm = 0
			for taxID in all_taxIDs:
				#note that taxID's might repeat in rank_taxIDs1 due to missing organism, so sum up all of these before taking the difference
				sum1=0
				for i in range(0,len(rank_taxIDs1)):
					if rank_taxIDs1[i]==taxID:
						sum1 = sum1 + rank_freqs1[i]
				sum2=0;
				for i in range(0,len(rank_taxIDs2)):
					if rank_taxIDs2[i]==taxID:
						sum2 = sum2 + rank_freqs2[i]
				L1norm = L1norm + abs(sum1-sum2)
		else:
			TP=-1
			FP=-1
			FN=-1
			Prec=-1
			Sens=-1
			L1norm=-1
		L1norms.append(L1norm)
		TPs.append(TP)
		FPs.append(FP)
		FNs.append(FN)
		Precs.append(Prec)
		Senss.append(Sens)
	
	#EMD code
	res = subprocess.check_output(python_loc + " " + EMD_loc + " -g " + input_file1 + " -r " + input_file2, shell=True)
	
	
	fid.write("TP \t" + '|'.join([str(TPs[i]) for i in range(0,len(rank_names))]) + "\n")
	fid.write("FP \t" + '|'.join([str(FPs[i]) for i in range(0,len(rank_names))]) + "\n")
	fid.write("FN \t" + '|'.join([str(FNs[i]) for i in range(0,len(rank_names))]) + "\n")
	fid.write("Prec \t" + '|'.join([str(Precs[i]) for i in range(0,len(rank_names))]) + "\n")
	fid.write("Sens \t" + '|'.join([str(Senss[i]) for i in range(0,len(rank_names))]) + "\n")
	fid.write("L1norm \t" + '|'.join([str(L1norms[i]) for i in range(0,len(rank_names))]) + "\n")
	fid.write("Unifrac \t" + res)
	fid.close()
	

if __name__ == "__main__":
	main(sys.argv[1:])

