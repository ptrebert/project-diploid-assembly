#!/usr/bin/env python

import os
import sys
import argparse
import glob


COLLECT_PATHS = {
    12: {
        'RUN_CONFIG': [
            '*.json',
            False,
            ''
        ],
        'READ_STATS_SUMMARY': [
            'input/fastq/*.stats',
            False,
            'statistics/input_reads'
        ],
        'CLUST_CONTIG_REF_BED': [
            'output/alignments/contigs_to_reference/reference_assembly/clustered/*/*.bed',
            True,
            'contig_coverage'
        ],
        'HAP_CONTIG_REF_BED': [
            'output/alignments/contigs_to_reference/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/polishing/*/haploid_assembly/*.bed',
            True,
            'contig_coverage'
        ],
        'HAP_CLUST_CONTIG_REF_BED': [
            'output/alignments/contigs_to_reference/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/polishing/*/clustering/*.fixed.bed',
            True,
            'contig_coverage/do_not_use'
        ],
        'COV_TRACKS': [
            'output/cov_tracks/hap_reads/longshot_QUAL10_GQ100/*/*/*/*.bigWig',
            True,
            'hap_read_coverage'
        ],
        'HAP_READS': [
            'output/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/draft/haploid_fastq/*h?-un*.fastq.gz',
            True,
            'haploid_reads'
        ],
        'HAP_READS_STATS_SUMMARY': [
            'output/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/draft/haploid_fastq/*h?-un*.stats',
            True,
            'statistics/haploid_reads/cluster'
        ],
        'HAP_READS_STATS_JOINT_SUMMARY': [
            'output/diploid_assembly/strandseq_joint/longshot_QUAL10_GQ100/*/*/*/draft/haploid_fastq/*.stats',
            True,
            'statistics/haploid_reads'
        ],
        'HAPLOTAGS': [
            'output/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/draft/haplotags/*.tsv',
            True,
            'haplotags'
        ],
        'HAPLOID_ASSEMBLIES': [
            'output/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/polishing/*/haploid_assembly/*h?-un.arrow-p1.fasta*',
            True,
            'assemblies/phased'
        ],
        'HAPLOID_CLUSTERED_ASSEMBLIES': [
            'output/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/polishing/*/clustering/*.fasta',
            True,
            'assemblies/do_not_use'
        ],
        'QUAST_REPORT_NHR_ASSEMBLY_TXT': [
            'output/evaluation/quastlg_busco/GRCh38_GCA_p13-GRCh38_GENCODEv31_basic/reference_assembly/non-hap-res/*/report.txt',
            False,
            'reports/quast'
        ],
        'QUAST_REPORT_NHR_ASSEMBLY_PDF': [
            'output/evaluation/quastlg_busco/GRCh38_GCA_p13-GRCh38_GENCODEv31_basic/reference_assembly/non-hap-res/*/report.pdf',
            False,
            'reports/quast'
        ],
        'QUAST_REPORT_CLUST_ASSEMBLY_TXT': [
            'output/evaluation/quastlg_busco/GRCh38_GCA_p13-GRCh38_GENCODEv31_basic/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/polishing/*/haploid_assembly/*/report.txt',
            True,
            'reports/quast'
        ],
        'QUAST_REPORT_CLUST_ASSEMBLY_PDF': [
            'output/evaluation/quastlg_busco/GRCh38_GCA_p13-GRCh38_GENCODEv31_basic/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/polishing/*/haploid_assembly/*/report.pdf',
            True,
            'reports/quast'
        ],
        'VARIANT_CALLS': [
            'output/integrative_phasing/longshot_QUAL10_GQ100/*/*/*/*.vcf.bgz*',
            True,
            'variant_calls'
        ],
        'PLOT_INPUT_READS': [
            'output/plotting/statistics/input_reads/*.stats.pdf',
            False,
            'plots/input_reads'
        ],
        'PLOT_SAARCLUST_DIAG_CLUST': [
            'output/plotting/saarclust_diagnostics/reference_assembly/clustered/*/*.pdf',
            True,
            'plots/saarclust_diagnostics'
        ],
        'PLOT_SAARCLUST_DIAG_HAP_CLUST': [
            'output/plotting/saarclust_diagnostics/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/polishing/*/clustering/*.pdf',
            True,
            'plots/saarclust_diagnostics/do_not_use'
        ],
        'NHR_ASSEMBLY': [
            'output/reference_assembly/non-hap-res/*.fasta*',
            False,
            'assemblies/unphased'
        ],
        'CLUST_NHR_ASSEMBLY': [
            'output/reference_assembly/clustered/*/*.fasta*',
            True,
            'assemblies/unphased'
        ],
        'CLUST_CONTIG_REF_STATS': [
            'output/statistics/contigs_to_ref_aln/reference_assembly/clustered/*/*.stats',
            True,
            'statistics/contig_alignments'
        ],
        'HAP_CONTIG_REF_STATS': [
            'output/statistics/contigs_to_ref_aln/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/polishing/*/haploid_assembly',
            True,
            'statistics/contig_alignments'
        ],
        'PHASING_STATS': [
            'output/statistics/phasing/longshot_QUAL10_GQ100/NA20847_hgsvc_pbsq2-clr_1000_scV12-hhu26/NA20847_hgsvc_pbsq2-clr_1000/NA20847_hgsvc_ilnxs-80pe_sseq/*wh-phased*',
            True,
            'statistics/phasing'
        ],
        'READ_STATS_DUMP': [
            'output/statistics/stat_dumps/*.pck',
            False,
            'statistics/input_reads/dumps'
        ],
        'HAP_READS_STATS_DUMPS': [
            'output/statistics/stat_dumps/diploid_assembly/strandseq_split/longshot_QUAL10_GQ100/*/*/*/draft/haploid_fastq/*.pck',
            True,
            'statistics/haploid_reads/cluster/dumps'
        ],
        'HAP_READS_STATS_JOINT_DUMPS': [
            'output/statistics/stat_dumps/diploid_assembly/strandseq_joint/longshot_QUAL10_GQ100/*/*/*/draft/haploid_fastq/*.pck',
            True,
            'statistics/haploid_reads/dumps'
        ],
        'HAPLOTAGGING_STATS': [
            'output/statistics/tag_split/longshot_QUAL10_GQ100/*/*/*/*.tsv',
            True,
            'statistics/haplotagging'
        ],
        'VARIANT_CALL_STATS': [
            'output/statistics/variant_calls/longshot/NA20847_hgsvc_pbsq2-clr_1000_scV12-hhu26/NA20847_hgsvc_ilnxs-80pe_sseq/*QUAL10.GQ100*.stats',
            True,
            'statistics/variant_calls'
        ]
    }
}


def parse_command_line():
    """
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--pipeline-work-dir',
        '-wd',
        type=str,
        dest='work_dir',
        default=os.getcwd(),
        help='Pipeline (Snkemake) working directory',
    )
    parser.add_argument(
        '--sample',
        '-s',
        type=str,
        dest='sample',
        required=True,
        help='Sample identifier, e.g., HG00733'
    )
    parser.add_argument(
        '--parameter-version',
        '-pv',
        type=int,
        dest='param_version',
        required=True,
        help='Parameter version, e.g., 12'
    )
    parser.add_argument(
        '--force-copy',
        '-fc',
        action='store_true',
        default=False,
        dest='force_copy',
        help='Copy files instead of symlinking'
    )
    parser.add_argument(
        '--dry-run',
        '-n',
        action='store_true',
        default=False,
        dest='dry_run',
        help='Dry run, print info but do not do anything'
    )
    parser.add_argument(
        '--dest-folder',
        '-d',
        type=str,
        dest='dest_folder',
        required=True,
        help='Path to folder for linking/copying. Will be created.'
    )
    parser.add_argument(
        '--file-report',
        '-f',
        type=str,
        dest='file_report',
        default=os.path.join(os.getcwd(), '{SAMPLE}.file_report.tsv'),
        help='Path to dump file report, defaults to current work dir.'
    )
    args = parser.parse_args()
    return args


def adapt_quast_report_name(file_path, output_path):
    """
    :param file_path:
    :param output_path:
    :return:
    """
    _, last_subdir = os.path.split(os.path.basename(file_path))
    out_name = last_subdir + '.quast-' + os.path.basename(file_path)
    return os.path.join(output_path, out_name)


def collect_result_files(top_source_path, sample, version, file_patterns, top_dest_path, dry_run):
    """
    :param top_source_path:
    :param sample:
    :param version:
    :param file_patterns:
    :param top_dest_path:
    :param dry_run::
    :return:
    """
    file_pairs = []
    for pattern_name, pattern_info in file_patterns.items():
        sub_pattern, check_version, out_folder = pattern_info
        out_folder = os.path.join(top_dest_path, out_folder)
        if not dry_run:
            os.makedirs(out_folder, exist_ok=True)
        pattern = os.path.join(top_source_path, sub_pattern)
        all_files = glob.glob(pattern)
        if not all_files:
            print('{}: 0 matches'.format(pattern_name))
            continue
        # sample ID must be somewhere
        all_files = [f for f in all_files if sample in f]
        if check_version:
            # version info must be somewhere
            all_files = [f for f in all_files if version in f]
        if not all_files:
            print('{}: 0 matches'.format(pattern_name))
            continue
        print('{}: {} matches'.format(pattern_name, len(all_files)))
        for src_file in all_files:
            if os.path.islink(src_file):
                src_file = os.readlink(src_file)
            if 'QUAST' in pattern_name:
                dest_file = adapt_quast_report_name(src_file, out_folder)
            else:
                src_name = os.path.basename(src_file)
                if out_folder.endswith('do_not_use'):
                    src_name = 'DEV-ONLY_'
                dest_file = os.path.join(out_folder, src_name)
            file_pairs.append((pattern_name, src_file, dest_file))
    print('Processing {} file pairs'.format(len(file_pairs)))
    return sorted(file_pairs)


def link_or_copy(file_pairs):
    """
    :param file_pairs:
    :return:
    """
    for _


def main():

    args = parse_command_line()
    version_string = 'scV{}'.format(args.param_version)
    glob_patterns = COLLECT_PATHS[args.param_version]
    print('Sample: {}'.format(args.sample))
    print('Pipeline parameter version: {} / {}'.format(args.param_version, version_string))
    pipeline_wd = os.path.abspath(args.work_dir)
    file_pairs = collect_result_files(
        pipeline_wd,
        args.sample,
        version_string,
        glob_patterns,
        args.dest_folder,
        args.dry_run
    )

    return 0

if __name__ == '__main__':
    sys.exit(main())
