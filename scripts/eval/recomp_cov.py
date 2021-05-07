#!/usr/bin/env python3

import os
import subprocess
import multiprocessing as mp

# input_bed = '/gpfs/project/ebertp/projects/rfdga/production/EVAL/run_folder/output/evaluation/hap_read_coverage'

# input_regions = '/gpfs/project/ebertp/projects/rfdga/production/EVAL/run_folder/output/evaluation/hap_read_coverage/temp'

# bigwig_files = sorted([os.path.join(input_bed, f) for f in os.listdir(input_bed) if f.endswith('bigWig')])

# region_files = sorted([os.path.join(input_regions, f) for f in os.listdir(input_regions) if f.endswith('loh.bed')])

# out_path = '/gpfs/project/ebertp/projects/rfdga/production/EVAL/run_folder/output/evaluation/hap_read_coverage/temp/avg_cov'

# for reg in region_files:
#     sample = os.path.basename(reg).rsplit('_', 1)[0]

#     for bigwig in bigwig_files:        
#         if sample not in bigwig:
#             continue
#         if 'h1-un' in bigwig or 'h2-un' in bigwig:
#             continue

#         bigwig_file = os.path.basename(bigwig).rsplit('.', 1)[0]

#         out_file = bed_file + 'cov_loh.tsv'
#         output = os.path.join(out_path, out_file)
        
#         shell_cmd = 'bigWigAverageOverBed {} {} {}'.format(bigwig, bed, output)

#         out = subprocess.check_output(shell_cmd, shell=True)
#         print('Done ', out_file)

###
# Below: compute assembly coverage in various GRCh38 regions
###

# input_bed = '/home/local/work/data/hgsvc/contig_aln_bed'

# input_regions = '/home/local/work/code/github/project-diploid-assembly/annotation/grch38'

# bed_files = sorted([os.path.join(input_bed, f) for f in os.listdir(input_bed) if f.endswith('.q60mrg.bed')])

# region_files = sorted([os.path.join(input_regions, f) for f in os.listdir(input_regions) if f.endswith('.regions')])

# out_path = '/home/local/work/data/hgsvc/aln_summary/regions'

# for reg in region_files:
#     reg_file = os.path.basename(reg).rsplit('.', 1)[0]

#     for bed in bed_files:
#         bed_file = os.path.basename(bed).rsplit('.', 1)[0]

#         out_file = bed_file + '_cov-in_' + reg_file + '.tsv'
#         output = os.path.join(out_path, out_file)
        
#         shell_cmd = 'bedtools coverage -a {} -b {} > {}'.format(reg, bed, output)

#         out = subprocess.check_output(shell_cmd, shell=True)
#         print('Done ', out_file)


###
# Below: compute assembly coverage in ROI (HLA, 3q29)
# use unfiltered/unmerged alignments to get count of 
# individual contigs
###

inv_comp_t2t = False

input_bed = '/home/local/work/data/hgsvc/contig_aln_bed/unfiltered/*.bed'

if inv_comp_t2t:
    # single analysis for T2T hap. assm. alignments
    input_bed = '/home/local/work/data/hgsvc/contig_aln_bed_t2t/*.bed'

# below: version for submission, MHC/HLA and 3q29
#input_regions = '/home/local/work/data/hgsvc/roi/roi_HLA_3q29.bed'
# below: version for revision, includes HLAn
#input_regions = '/home/local/work/data/hgsvc/roi/roi_HLA_3q29_HLAn.bed'
# below: version post publication, for inversion companion, includes 1p36.13
input_regions = '/home/local/work/data/hgsvc/roi/roi_HLA_3q29_HLAn_1p3613.bed'

if inv_comp_t2t:
    input_regions = '/home/local/work/data/hgsvc/roi/roi-t2t_1p3613.bed'

#output = '/home/local/work/data/hgsvc/roi/MAPQ-ALL_roi_assm_overlaps.tsv'
#output = '/home/local/work/data/hgsvc/roi/MAPQ-ALL_roi_HLAn_assm_overlaps.tsv'
output = '/home/local/work/data/hgsvc/roi/MAPQ-ALL_roi_HLAn_1p3613_assm_overlaps.tsv'

if inv_comp_t2t:
    output = '/home/local/work/data/hgsvc/roi/MAPQ-ALL_roi-t2t_1p3613_assm_overlaps.tsv'

shell_cmd = 'bedtools intersect -filenames -wao -a {} -b {} > {}'.format(input_regions, input_bed, output)

out = subprocess.check_output(shell_cmd, shell=True)
print('Done ', output)


def execute_system_calls(calls):

    for c in calls:
        if c is None:
            continue
        print('exec ', c)
        out = subprocess.check_output(c, shell=True)
    return 0

# ###
# # Below: make dot plots for chr3q29
# ###

# path = '/home/local/work/data/hgsvc/roi/dotplots_chr3'
# ref_fasta = os.path.join(path, 'GRCh38_HGSVC2_noalt.chr3.fasta')

# scaffolds = [f for f in os.listdir(path) if 'GRCh38' not in f and f.endswith('.fasta')]

# exec_commands = []

# samples = [
#     ('HG00096', 'CLR', 'H2'),
#     ('HG00512', 'CLR', 'H2'),
#     ('HG00514', 'CLR', 'H2'),
#     ('HG00733', 'CLR', 'H1'),
#     ('HG00733', 'CLR', 'H2'),
#     ('HG00864', 'CLR', 'H2'),
#     ('HG02492', 'CLR', 'H1'),
#     ('HG03065', 'CLR', 'H1'),
#     ('HG03065', 'CLR', 'H2'),
#     ('HG03371', 'CLR', 'H2'),
#     ('NA19650', 'CLR', 'H2'),
#     ('NA20509', 'CLR', 'H1')
# ]

# for scf in scaffolds:
#     sample, _, platform, _ = scf.split('_', 3)
#     tech = platform.split('-')[1].upper()
#     if 'h1-un' in scf:
#         hap = 'H1'
#     else:
#         hap = 'H2'
#     if (sample, tech, hap) not in samples:
#         continue
#     query_fasta = os.path.join(path, scf)
#     min_length = 1000
#     out_mums = os.path.join(path, 'mums', '{}_{}_{}.{}.chr3.mums'.format(sample, tech, hap, min_length))

#     if not os.path.isfile(out_mums):
#         mum_call = 'mummer -mum -n -c -b -l {} -F {} {} > {}'.format(min_length, ref_fasta, query_fasta, out_mums)
#     else:
#         mum_call = None
#         #out = subprocess.check_output(mum_call, shell=True)

#     plot_prefix = os.path.join(path, 'plots', '{}_{}_{}.{}.chr3'.format(sample, tech, hap, min_length))
#     plot_call = 'mummerplot -s medium -p {} -t png -R {} -Q {} {}'.format(plot_prefix, ref_fasta, query_fasta, out_mums)
#     #out = subprocess.check_output(plot_call, shell=True)

#     plot_script = plot_prefix + '.gp'
#     gnuplot_call = 'gnuplot {}'.format(plot_script)
#     #out = subprocess.check_output(gnuplot_call, shell=True)
#     exec_commands.append((mum_call, plot_call, gnuplot_call))

# ###
# # Below: align Bionano segments against phased assembly contigs to confirm placement
# ###


# path = '/home/local/work/data/hgsvc/roi/dotplots_chr3/chr3_ctg_fasta'
# ref_fasta_path = '/home/local/work/data/hgsvc/roi/dotplots_chr3/chr3_segments'
# out_path = '/home/local/work/data/hgsvc/roi/dotplots_chr3/chr3_alignments/contigs'
# #ref_fastas = sorted([os.path.join(ref_fasta_path, f) for f in os.listdir(ref_fasta_path)])
# ref_fasta = os.path.join(ref_fasta_path, 'hg38_chr3q29_all_segments.fasta')

# contigs = [f for f in os.listdir(path) if 'GRCh38' not in f and f.endswith('ctg3q29.fasta')]

# minimap_call = 'minimap2 -t 1 --secondary=yes --eqx -Y -ax asm20 -R "{readgroup}" {reference} {query} | samtools sort | samtools view -b -F 260 /dev/stdin > {output}'
# bedtools_call = 'bedtools bamtobed -i {input} > {output}'

# exec_commands = []

# for ctg in contigs:
#     sample, _, platform, _ = ctg.split('_', 3)
#     tech = platform.split('-')[1].upper()
#     if 'h1-un' in ctg:
#         hap = 'H1'
#     else:
#         hap = 'H2'

#     contig_fasta = os.path.join(path, ctg)
#     out_bam = os.path.join(out_path, '{}_{}_{}.all-blocks.bam'.format(sample, hap, tech))
#     out_bed = os.path.join(out_path, '{}_{}_{}.all-blocks.bed'.format(sample, hap, tech))
#     run_minimap = None
#     run_bedtools = None
#     if not os.path.isfile(out_bam):
#         readgroup = '@RG\\tID:{}\\tSM:{}'.format('_'.join([sample, hap, tech]), sample)
#         run_minimap = minimap_call.format(**{'readgroup': readgroup, 'reference': contig_fasta, 'query': ref_fasta, 'output': out_bam})
#     if not os.path.isfile(out_bed):
#         run_bedtools = bedtools_call.format(**{'input': out_bam, 'output': out_bed})
#     exec_commands.append((run_minimap, run_bedtools))

# paths = [
#     '/home/local/work/data/hgsvc/contig_aln_bed/unfiltered/grt00',
#     '/home/local/work/data/hgsvc/contig_aln_bed/unfiltered/grt20'
# ]
# paths = [
#     '/home/local/work/data/hgsvc/contig_aln_bed_t2t/grt00',
#     '/home/local/work/data/hgsvc/contig_aln_bed_t2t/grt20'
# ]
# genome_file = '/home/local/work/data/hgsvc/noalt_ref/hg38.no_alt.fa.gz.fai'
# genome_file = '/home/local/work/data/hgsvc/t2tv1/T2Tv1_38p13Y_chm13.fasta.fai'

# exec_commands = []

# for p in paths:
#     bed_files = [os.path.join(p, b) for b in os.listdir(p) if b.endswith('.bed') and 'complement' not in b]
#     out_files = [b.replace('.bed', '.complement.bed') for b in bed_files]
#     for b, o in zip(bed_files, out_files):
#         syscall = 'bedtools complement -L -i {} -g {} > {}'.format(b, genome_file, o)
#         exec_commands.append((syscall, ))


# with mp.Pool(4) as pool:
#     resit = pool.imap_unordered(execute_system_calls, exec_commands)
#     print('processing')
#     for res in resit:
#         print('done')
