#!/usr/bin/env python3

"""
This script realizes the comparative assembly analysis
as described in the manuscript. For questions about
this analysis or this script in particular, please
contact the contributors
Jana Ebler or Tobias Marschall
(see https://doi.org/10.1101/855049 )
"""

import sys
from collections import defaultdict, namedtuple
from pysam import VariantFile
from whatshap.args import HelpfulArgumentParser as ArgumentParser


def check_mendel(mother_gt, father_gt, child_gt):
	return ((child_gt[0] in mother_gt) and (child_gt[1] in father_gt)) or \
		((child_gt[0] in father_gt) and (child_gt[1] in mother_gt))


def main():
	parser = ArgumentParser(prog='compare-vcfs.py', description=__doc__)
	parser.add_argument('--parent-names', default=None,
		help='Names of parent samples (comma separated) for Mendelian check')
	parser.add_argument('vcf', metavar='VCF', help='VCF files to compare')
	parser.add_argument('samples', nargs='+', metavar='SAMPLES', help='Samples to compare')
	args = parser.parse_args()

	vcf = VariantFile(args.vcf)
	
	mother, father = None, None
	has_parents = False
	if args.parent_names:
		mother, father = args.parent_names.split(',')
		has_parents = True
		mendelian_errors = defaultdict(int)
	
	#print(dir(vcf))
	
	n = 0
	all_concordant_count = 0
	fully_discordant_count = 0
	minority_errors = defaultdict(int)
	CallReport = namedtuple('CallReport', ['mendel_ok','minority_gt','all_concordant', 'fully_discordant'])
	
	call_reports = dict()
	for sample in args.samples:
		call_reports[sample] = defaultdict(int)

	for record in vcf:
		n += 1
		#print(dir(record))
		#print(list(record.samples))

		# Sort genotypes (i.e. disregard phase)
		gts = [tuple(sorted(record.samples[sample]['GT'])) for sample in args.samples]

		mendel_ok_list = [None] * len(gts)
		minority_list = [None] * len(gts)

		if has_parents:
			mother_gt = record.samples[mother]['GT']
			father_gt = record.samples[father]['GT']
		for i, sample in enumerate(args.samples):
			call = record.samples[sample]
			gt = call['GT']
			mendel_ok = check_mendel(mother_gt,father_gt,gt)
			mendel_ok_list[i] = mendel_ok
			if not mendel_ok:
				mendelian_errors[sample] += 1

		all_concordant = None
		fully_discordant = None
		#print(gts)
		if len(gts) > 1:
			gt_counts = defaultdict(int)
			for gt in gts:
				gt_counts[gt] += 1
			majority_count, majority_gt = sorted( ((count,gt) for gt,count in gt_counts.items()), reverse=True)[0]
			# if there is no majority, use the first genotype
			if majority_count <= len(gts)/2:
				majority_gt = gts[0]
			if has_parents:
				if not check_mendel(mother_gt,father_gt,majority_gt):
					mendelian_errors['MAJORITY'] += 1
				
			all_concordant = len(gt_counts) == 1
			fully_discordant = len(gt_counts) > 2
			#print(gts, gt_counts)
			if all_concordant:
				all_concordant_count += 1
				minority_list = [False] * len(gts)
			elif fully_discordant:
				fully_discordant_count += 1
			else:
				for i, (sample, gt) in enumerate(zip(args.samples, gts)):
					minority_list[i] = gt != majority_gt 
					if gt != majority_gt:
						minority_errors[sample] += 1

		for i, (sample, mendel_ok, minority_gt) in enumerate(zip(args.samples, mendel_ok_list, minority_list)):
			report = CallReport(mendel_ok, minority_gt, all_concordant, fully_discordant)
			call_reports[sample][report] += 1
	
	print('Processed {} variant sites'.format(n))
	print('Concordant across all samples: {}, {:.2f}%'.format(all_concordant_count, all_concordant_count/n*100))
	print('Fully discordant (more than two different genotypes observed): {}, {:.2f}%'.format(fully_discordant_count, fully_discordant_count/n*100))
	print('Minority genotypes (likely errors):')
	for sample in args.samples:
		print('\t{}: {} {:.2f}%'.format(sample, minority_errors[sample], minority_errors[sample]/n*100))
	if has_parents:
		print('Mendelian errors:')
		for sample in mendelian_errors.keys():
			print('\t{}: {} {:.2f}%'.format(sample, mendelian_errors[sample], mendelian_errors[sample]/n*100))
	print('Call reports:')
	for sample in args.samples:
		print('\t', sample)
		for count, report in sorted(((c, r) for r, c in call_reports[sample].items()), reverse=True):
			print('\t\t', report, count)
		
		#print(sample, call_reports[sample])
		#break


if __name__ == '__main__':
	main()
