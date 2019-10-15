
include: 'smk_include/handle_data_download.smk'
include: 'smk_include/prepare_custom_references.smk'
include: 'smk_include/results_child.smk'
include: 'smk_include/results_parents.smk'

localrules: master

# Global wildcard constraints for consistent naming
wildcard_constraints:
    # Reference genome or assembly
    reference = '[A-Za-z0-9_\-]+',
    # STS = STrand-Seq reads used for clustering and integrative phasing
    sts_reads = '[A-Za-z0-9_\-]+',
    # VC = variant calling, read set used for variant calling
    vc_reads = '[A-Za-z0-9_\-]+',
    # HAP = haplotype, read set that is tagged and split by haplotype
    hap_reads = '[A-Za-z0-9_\-]+',
    # POL = polishing, read set used for polishing (aligned against contigs to be polished)
    pol_reads = '[A-Za-z0-9_\-]',
    # allowed variant callers
    var_caller = '(freebayes|longshot|deepvar)',
    # allowed assembly tools
    assembler = '(wtdbg|canu|flye|pereg)',
    # GQ / DP = genotype quality / depth at position
    gq = '[0-9]+',
    dp = '[0-9]+',
    # haplotype identifier
    hap = '[h12un\-]+',
    # sequence = chromosome, contig, cluster etc.
    sequence = '[A-Za-z0-9]+'


rule master:
    input:
        # this triggers a checkpoint
        # for downloading the strand-seq data
        expand('input/fastq/strand-seq/{individual}_{bioproject}/requests',
                individual=['HG00733', 'HG00732', 'HG00731'],
                bioproject=['PRJEB12849']),
        rules.master_results_child.input,
        rules.master_results_parents.input

    message: 'Executing ALL'


def make_log_useful(log_path, status):

    my_env = dict(os.environ)
    with open(log_path, 'a') as logfile:
        _ = logfile.write('\n===[{}]===\n'.format(status))
        _ = logfile.write('Host: {}\n'.format(my_env.get('HOSTNAME', 'N/A')))
        _ = logfile.write('Display: {}\n'.format(my_env.get('DISPLAY', 'N/A')))
        _ = logfile.write('Shell: {}\n'.format(my_env.get('SHELL', 'N/A')))
        _ = logfile.write('Terminal: {}\n'.format(my_env.get('TERM', 'N/A')))
        _ = logfile.write('Screen: {}\n'.format(my_env.get('STY', 'N/A')))
        _ = logfile.write('Conda ENV: {}\n'.format(my_env.get('CONDA_DEFAULT_ENV', 'N/A')))
        _ = logfile.write('\n')
    return


onsuccess:
    make_log_useful(log, 'SUCCESS')
    import socket
    host = socket.gethostname()
    if config['notify'] and 'bibigrid' not in host:
        shell('mail -s "[Snakemake] DGA - SUCCESS" {} < {{log}}'.format(config['notify_email']))


onerror:
    make_log_useful(log, 'ERROR')
    import socket
    host = socket.gethostname()
    if config['notify'] and 'bibigrid' not in host:
        shell('mail -s "[Snakemake] DGA - ERRROR" {} < {{log}}'.format(config['notify_email']))
