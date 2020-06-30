#!/usr/bin/env python3

import os
import sys
import argparse
import glob
import shutil
import errno
import collections


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
            'output/alignments/contigs_to_reference/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/haploid_assembly/*.bed',
            True,
            'contig_coverage'
        ],
        'HAP_CLUST_CONTIG_REF_BED': [
            'output/alignments/contigs_to_reference/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/clustering/*.fixed.bed',
            True,
            'contig_coverage/do_not_use'
        ],
        'COV_TRACKS': [
            'output/cov_tracks/hap_reads/*/*/*/*/*.bigWig',
            True,
            'hap_read_coverage'
        ],
        'HAP_READS': [
            'output/diploid_assembly/strandseq_split/*/*/*/*/draft/haploid_fastq/*h?-un*.fastq.gz',
            True,
            'haploid_reads'
        ],
        'HAP_READS_STATS_SUMMARY': [
            'output/diploid_assembly/strandseq_split/*/*/*/*/draft/haploid_fastq/*h?-un*.stats',
            True,
            'statistics/haploid_reads/cluster'
        ],
        'HAP_READS_STATS_JOINT_SUMMARY': [
            'output/diploid_assembly/strandseq_joint/*/*/*/*/draft/haploid_fastq/*.stats',
            True,
            'statistics/haploid_reads'
        ],
        'HAPLOTAGS': [
            'output/diploid_assembly/strandseq_split/*/*/*/*/draft/haplotags/*.tsv',
            True,
            'haplotags'
        ],
        'HAPLOID_ASSEMBLIES': [
            [
                'output/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/haploid_assembly/*h?-un.arrow-p1.fasta*',
                'output/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/haploid_assembly/*h?-un.racon-p2.fasta*'
            ],
            True,
            'assemblies/phased'
        ],
        'HAPLOID_CLUSTERED_ASSEMBLIES': [
            'output/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/clustering/*.fasta',
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
            'output/evaluation/quastlg_busco/GRCh38_GCA_p13-GRCh38_GENCODEv31_basic/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/haploid_assembly/*/report.txt',
            True,
            'reports/quast'
        ],
        'QUAST_REPORT_CLUST_ASSEMBLY_PDF': [
            'output/evaluation/quastlg_busco/GRCh38_GCA_p13-GRCh38_GENCODEv31_basic/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/haploid_assembly/*/report.pdf',
            True,
            'reports/quast'
        ],
        'VARIANT_CALLS': [
            'output/integrative_phasing/*/*/*/*/*.vcf.bgz*',
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
            'output/plotting/saarclust_diagnostics/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/clustering/*.pdf',
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
            'output/statistics/contigs_to_ref_aln/diploid_assembly/strandseq_split/*/*/*/*/polishing/*/haploid_assembly/*.stats',
            True,
            'statistics/contig_alignments'
        ],
        'PHASING_STATS': [
            'output/statistics/phasing/*/*/*/*/*wh-phased*',
            True,
            'statistics/phasing'
        ],
        'READ_STATS_DUMP': [
            'output/statistics/stat_dumps/*.pck',
            False,
            'statistics/input_reads/dumps'
        ],
        'HAP_READS_STATS_DUMPS': [
            'output/statistics/stat_dumps/diploid_assembly/strandseq_split/*/*/*/*/draft/haploid_fastq/*.pck',
            True,
            'statistics/haploid_reads/cluster/dumps'
        ],
        'HAP_READS_STATS_JOINT_DUMPS': [
            'output/statistics/stat_dumps/diploid_assembly/strandseq_joint/*/*/*/*/draft/haploid_fastq/*.pck',
            True,
            'statistics/haploid_reads/dumps'
        ],
        'HAPLOTAGGING_STATS': [
            'output/statistics/tag_split/*/*/*/*/*.tsv',
            True,
            'statistics/haplotagging'
        ],
        'VARIANT_CALL_STATS': [
            'output/statistics/variant_calls/*/*/*/*QUAL10.GQ100*.stats',
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
        '--ignore-existing-dest',
        '-igex',
        action='store_true',
        default=False,
        dest='ignore_existing',
        help='Ignore files already existing in "dest-folder"'
    )
    parser.add_argument(
        '--ignore-missing-files',
        '-igms',
        action='store_true',
        default=False,
        dest='ignore_missing',
        help='Ignore missing result files.'
    )
    parser.add_argument(
        '--show-glob-patterns',
        '-globs',
        action='store_true',
        default=False,
        dest='show_globs',
        help='Just print glob patterns and exist'
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


def adapt_quast_report_name(file_path, output_path, version_prefix):
    """
    :param file_path:
    :param output_path:
    :return:
    """
    _, last_subdir = os.path.split(os.path.dirname(file_path))
    out_name = last_subdir + '.quast-' + os.path.basename(file_path)
    return os.path.join(output_path, version_prefix + out_name)


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
    version_prefix = 'v{}_'.format(version)
    version_string = 'scV{}'.format(version)
    count_matches = dict()
    file_pairs = []
    for pattern_name, pattern_info in file_patterns.items():
        sub_patterns, check_version, out_folder = pattern_info
        out_folder = os.path.join(top_dest_path, out_folder)
        if not dry_run:
            os.makedirs(out_folder, exist_ok=True)
        if isinstance(sub_patterns, str):
            sub_patterns = [sub_patterns]
        file_pattern_count = 0
        for sub_pattern in sub_patterns:
            pattern = os.path.join(top_source_path, sub_pattern)
            all_files = glob.glob(pattern)
            if not all_files:
                continue
            # sample ID must be somewhere
            all_files = [f for f in all_files if sample in f]
            if check_version:
                # version info must be somewhere
                all_files = [f for f in all_files if version_string in f]
            if not all_files:
                continue
            file_pattern_count += len(all_files)
            for src_file in all_files:
                src_name = None
                if os.path.islink(src_file):
                    src_name = os.path.basename(src_file)
                    src_file = os.readlink(src_file)
                if 'QUAST' in pattern_name:
                    dest_file = adapt_quast_report_name(src_file, out_folder, version_prefix)
                else:
                    if src_name is None:
                        src_name = os.path.basename(src_file)
                    if out_folder.endswith('do_not_use'):
                        src_name = 'DEV-ONLY_'.format(version_prefix) + src_name
                    dest_file = os.path.join(out_folder, version_prefix + src_name)
                assert os.path.isfile(src_file), 'Not a valid file path: {}'.format(src_file)
                file_pairs.append((pattern_name, src_file, dest_file))
        count_matches[(pattern_name, ' // '.join(sub_patterns))] = file_pattern_count
    print('Processing {} files'.format(len(file_pairs)))
    return sorted(file_pairs), count_matches


def link_or_copy(file_pairs, dry_run, force_copy, ignore_existing):
    """
    :param file_pairs:
    :param dry_run:
    :param force_copy:
    :param ignore_existing:
    :return:
    """
    for num, (_, src, dest) in enumerate(file_pairs, start=1):
        if dry_run:
            print('{}. - - - - - - - -'.format(num))
            print('SRC: {}'.format(src))
            print('DEST: {}'.format(dest))
            print('- - - - - - - - - -')
            continue
        try:
            if force_copy:
                shutil.copy(src, dest)
            else:
                os.symlink(src, dest)
        except OSError as e:
            if e.errno == errno.EEXIST:
                if ignore_existing:
                    print('WARN - destination exists: {}'.format(dest))
                    pass
                else:
                    raise e
    return


def main():

    args = parse_command_line()
    glob_patterns = COLLECT_PATHS[args.param_version]
    print('Sample: {}'.format(args.sample))
    print('Pipeline parameter version: {}'.format(args.param_version))

    if args.show_globs:
        print('Pipeline_folder\tDestination_folder\n')
        for n, _, p in sorted(glob_patterns.values()):
            print('{}\t{}'.format(n, p))
        return 0

    pipeline_wd = os.path.abspath(args.work_dir)
    file_pairs, count_matches = collect_result_files(
        pipeline_wd,
        args.sample,
        args.param_version,
        glob_patterns,
        args.dest_folder,
        args.dry_run
    )
    if not args.ignore_missing:
        missing_results = [k for k, v in count_matches.items() if v == 0]
        if missing_results:
            raise ValueError('Missing results:\n{}'.format('\n'.join(['{}\t{}'.format(n, g) for n, g in missing_results])))


    link_or_copy(file_pairs, args.dry_run, args.force_copy, args.ignore_existing)

    if not args.dry_run:
        file_report = args.file_report.format(**{'SAMPLE': args.sample})
        file_report = os.path.abspath(file_report)
        os.makedirs(os.path.dirname(file_report), exist_ok=True)

        with open(file_report, 'w') as dump:
            _ = dump.write('#RESULT_TYPE\tSRC_FILE\tDEST_FILE\n')
            _ = dump.write('\n'.join(['\t'.join(record) for record in file_pairs]))

    return 0


if __name__ == '__main__':
    sys.exit(main())
