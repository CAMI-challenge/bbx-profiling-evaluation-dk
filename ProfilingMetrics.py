#!/usr/bin/env python

import os
import sys
import getopt
from EMDUnifrac import emd_unifrac

# input_file1='../EMD/Test/test_truth_CAMI.txt'
# input_file2='../EMD/Test/test_reconstruction_CAMI.txt'
# output_file='test.txt'
# epsilon = 0
# normalize="n"


def main(argv):
	file_path_truth = None
	file_path_recon = None
	file_path_output = None
	epsilon = None
	normalize = True
	try:
		opts, args = getopt.getopt(
			argv, "g:r:o:e:n:h",
			["GroundTruth=", "Reconstruction=", "Output=", "Epsilon=", "Normalize="])
	except getopt.GetoptError:
		print 'Call using: python ProfilingMetrics.py -g <GroundTruth.profile> -r <Reconstruction.profile> -o <Output.txt> -e epsilon -n <normalize(y/n)'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Call using: python ProfilingMetrics.py -g <GroundTruth.profile> -r <Reconstruction.profile> -o <Output.txt> -e epsilon -n <normalize(y/n)'
			sys.exit(2)
		elif opt in ("-g", "--GroundTruth"):
			file_path_truth = arg
		elif opt in ("-r", "--Reconstruction"):
			file_path_recon = arg
		elif opt in ("-o", "--Output"):
			file_path_output = arg
		elif opt in ("-e", "--Epsilon"):
			epsilon = float(arg)
		elif opt in ("-n", "--Normalize"):
			normalize = arg == 'y'
	calc_metrics(file_path_truth, file_path_recon, file_path_output, epsilon=epsilon, normalize=normalize)


def read_taxonomy_file(file_path, epsilon=None):
	assert isinstance(file_path, basestring)
	assert epsilon is None or isinstance(epsilon, (float, int, long))
	tax_path = list()
	tax_ids = list()
	weights = dict()
	ranks = []
	with open(file_path, 'r') as read_handler:
		for line in read_handler:
			line = line.rstrip()
			if len(line) == 0:
				continue  # skip blank lines
			if line.lower().split(':')[0] == '@ranks':
				ranks = line.strip().split(':')[1].split('|')
				continue
			if line[0] in ['@', '#']:
				continue  # skip comment or header
			temp_split = line.split('\t')
			weight = float(temp_split[4])
			if epsilon is not None and weight < epsilon:
				# Ignore taxIDs below cutoff
				continue
			if weight == 0:
				# Ignore zero weighted taxIDs
				continue
			tax_path.append(temp_split[2])  # add the whole taxpath
			tax_ids.append(temp_split[0])  # just terminal tax ID
			weights[temp_split[0]] = weight  # the associated weight
	return tax_ids, tax_path, weights, ranks


def calc_metrics(file_path_truth, file_path_recon, file_path_output, epsilon=None, normalize=True):
	assert isinstance(file_path_truth, basestring), file_path_truth
	assert isinstance(file_path_recon, basestring), file_path_recon
	assert isinstance(file_path_output, basestring), file_path_output
	assert epsilon is None or isinstance(epsilon, (float, int, long)), epsilon
	assert isinstance(normalize, bool), normalize

	# read taxonomic profiles
	tax_ids1, tax_path1, weights1_dict, ranks1 = read_taxonomy_file(file_path_truth, epsilon)
	tax_ids2, tax_path2, weights2_dict, ranks2 = read_taxonomy_file(file_path_recon, epsilon)

	# Get the ranks, whichever is longest
	if len(ranks1) > len(ranks2):
		rank_names = ranks1
	else:
		rank_names = ranks2

	# open the output file
	fid = open(file_path_output, 'w')
	fid.write("@input\t " + os.path.basename(file_path_recon) + "\n")
	fid.write("@ground_truth\t " + os.path.basename(file_path_truth) + "\n")
	fid.write("@Taxonomic ranks\t " + '|'.join(rank_names) + "\n")

	list_tps = list()
	list_fps = list()
	list_fns = list()
	list_precs = list()
	list_senss = list()
	list_l1_norms = list()

	for rank_pos, rank in enumerate(rank_names):
		# select tax IDs and freqs of interest from input1
		# rank_pos = rank_names.index(rank)
		# get the tax ids at the rank of interest
		rank_tax_ids1 = list()
		rank_freqs1 = list()
		for i in range(0, len(tax_path1)):
			temp = tax_path1[i].split('|')
			if len(temp) == rank_pos + 1:
				rank_tax_ids1.append(temp[rank_pos])
				rank_freqs1.append(weights1_dict[temp[rank_pos]])
		# select tax IDs and freqs of interest from input2
		rank_tax_ids2 = list()
		rank_freqs2 = list()
		for i in range(0, len(tax_path2)):
			temp = tax_path2[i].split('|')
			if len(temp) == rank_pos + 1:
				rank_tax_ids2.append(temp[rank_pos])
				rank_freqs2.append(weights2_dict[temp[rank_pos]])

		# normalize if this was asked for
		if normalize:
			if sum(rank_freqs1) > 0:
				rank_freqs1 = [rank_freqs1[i] / sum(rank_freqs1) for i in range(0, len(rank_freqs1))]
			if sum(rank_freqs2) > 0:
				rank_freqs2 = [rank_freqs2[i] / sum(rank_freqs2) for i in range(0, len(rank_freqs2))]
		set_rank_tax_ids1 = set(rank_tax_ids1)
		set_rank_tax_ids2 = set(rank_tax_ids2)
		if sum(rank_freqs1) + sum(rank_freqs2) > 0:
			# TP/FP/etc
			true_positives = len(set_rank_tax_ids1 & set_rank_tax_ids2)  # in both
			false_positives = len(set_rank_tax_ids2.difference(set_rank_tax_ids1))  # in reconstruction not truth
			false_negatives = len(set_rank_tax_ids1.difference(set_rank_tax_ids2))  # in truth not reconstruction
			if true_positives + false_positives == 0:
				precision = -1
			else:
				precision = float(true_positives) / (true_positives + false_positives)
			if true_positives + false_negatives == 0:
				sensitivity = -1
			else:
				sensitivity = float(true_positives) / (true_positives + false_negatives)
			# L1 error
			all_tax_ids = set_rank_tax_ids1 | set_rank_tax_ids2
			l1_norm = 0
			for taxID in all_tax_ids:
				# note that taxID's might repeat in rank_taxIDs1 due to missing organism,
				# so sum up all of these before taking the difference
				sum1 = 0
				for i in range(0, len(rank_tax_ids1)):
					if rank_tax_ids1[i] == taxID:
						sum1 = sum1 + rank_freqs1[i]
				sum2 = 0
				for i in range(0, len(rank_tax_ids2)):
					if rank_tax_ids2[i] == taxID:
						sum2 = sum2 + rank_freqs2[i]
				l1_norm = l1_norm + abs(sum1 - sum2)
		else:
			true_positives = -1
			false_positives = -1
			false_negatives = -1
			precision = -1
			sensitivity = -1
			l1_norm = -1
		list_l1_norms.append(l1_norm)
		list_tps.append(true_positives)
		list_fps.append(false_positives)
		list_fns.append(false_negatives)
		list_precs.append(precision)
		list_senss.append(sensitivity)

	# EMD code
	res = emd_unifrac(file_path_truth, file_path_recon)
	# res = subprocess.check_output(python_loc + " " + EMD_loc + " -g " + input_file1 + " -r " + input_file2, shell=True)

	fid.write("TP \t" + '|'.join([str(list_tps[i]) for i in range(0, len(rank_names))]) + "\n")
	fid.write("FP \t" + '|'.join([str(list_fps[i]) for i in range(0, len(rank_names))]) + "\n")
	fid.write("FN \t" + '|'.join([str(list_fns[i]) for i in range(0, len(rank_names))]) + "\n")
	fid.write("Prec \t" + '|'.join([str(list_precs[i]) for i in range(0, len(rank_names))]) + "\n")
	fid.write("Sens \t" + '|'.join([str(list_senss[i]) for i in range(0, len(rank_names))]) + "\n")
	fid.write("L1norm \t" + '|'.join([str(list_l1_norms[i]) for i in range(0, len(rank_names))]) + "\n")
	fid.write("Unifrac \t" + str(res))
	fid.close()


if __name__ == "__main__":
	main(sys.argv[1:])
