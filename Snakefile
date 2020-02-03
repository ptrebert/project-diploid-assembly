
include: 'smk_include/module_includes.smk'

include: 'smk_include/results/run_afr_trios.smk'
include: 'smk_include/results/run_amr_trios.smk'
include: 'smk_include/results/run_eas_trios.smk'
include: 'smk_include/results/run_eur_trios.smk'
include: 'smk_include/results/run_sas_trios.smk'

localrules: master, setup_env

# Global wildcard constraints for consistent naming
wildcard_constraints:
    # For Strand-seq haploid DGA, possible to proceed in a per-haplotype/per-cluster fashion (or jointly)
    hap_assm_mode = '(joint|split)',
    # Sample should be used for individual + project + platform (+ x, but no file extension)
    sample = '[A-Za-z0-9_\-]+',
    # Individual such as HG00733
    individual = '[A-Z0-9]+',
    # Reference genome or assembly
    reference = '[A-Za-z0-9_\-]+',
    # STS = STrand-Seq reads used for clustering and integrative phasing
    sts_reads = '[A-Za-z0-9_\-]+',
    # VC = variant calling, read set used for variant calling
    vc_reads = '[A-Za-z0-9_\-]+',
    # HAP = haplotype, read set that is tagged and split by haplotype
    hap_reads = '[A-Za-z0-9_\-]+',
    # POL = polishing, read set used for polishing (aligned against contigs to be polished)
    pol_reads = '[A-Za-z0-9_\-]+',
    # allowed variant callers
    var_caller = '(freebayes|longshot|deepvar)',
    # allowed assembly tools
    assembler = '(wtdbg|canu|flye|pereg|shasta)',
    # polisher
    polisher = '(arrow|racon)',
    pol_pass ='(arrow|racon)\-p[12]',
    # GQ / DP = genotype quality / depth at position
    gq = '[0-9]+',
    dp = '[0-9]+',
    qual = '[0-9]+',
    # haplotype identifier
    hap = '[h12un\-]+',
    # sequence = chromosome, contig, cluster etc.
    sequence = '[A-Za-z0-9]+',
    # For files that may or may not be just a split (chrom, contig etc.)
    # of the complete file, can carry a qualifier after the "hap"
    split = '(^$|^\.[A-Za-z0-9]+)',
    # Some generic constraints to enforce a more consistent naming scheme
    folder_path = '[A-Za-z0-9\-_\/]+',  # note: "." is NOT allowed in a folder path
    file_name = '[A-Za-z0-9\-_\.]+',
    known_ref = 'GRCh3[78][A-Za-z0-9_]+',
    genemodel = 'GRCh38[A-Za-z0-9_]+'


rule master:
    input:
        rules.run_afr_trios.input,
        rules.run_amr_trios.input,
        rules.run_eur_trios.input,
        rules.run_eas_trios.input,
        rules.run_sas_trios.input
    message: 'Default run: processing all HGSVC samples'


def collect_all_configured_samples(wildcards):
    """
    :return:
    """
    configured_samples = []
    build_target_path = 'output/targets/{super_pop}_{population}_{family_id}/{individual}.fofn'

    for key, values in config.items():
        if not key.startswith('sample_description'):
            continue
        sample_desc = {
            'super_pop': values['super_population'],
            'population': values['population'],
            'family_id': values['family'],
            'individual': values['individual']
        }
        sample_targets = build_target_path.format(**sample_desc)
        configured_samples.append(sample_targets)
    configured_samples = sorted(set(configured_samples))
    return configured_samples


rule master_custom:
    input:
        collect_all_configured_samples
    message: 'Custom run: processing only samples in loaded configuration'


rule setup_env:
    input:
        rules.create_cluster_log_folders.output,
        rules.create_conda_environment_shell_tools.output,
        rules.create_conda_environment_pacbio_tools.output,
        rules.create_conda_environment_r_tools.output,
        rules.create_conda_environment_r_script.output,
        rules.create_conda_environment_bio_tools.output,
        rules.create_conda_environment_pyscript.output,
        rules.download_shasta_executable.output,
        rules.download_quast_busco_databases.output,
        rules.check_singularity_version.output,
        rules.install_rlib_saarclust.output,
        rules.install_rlib_breakpointr.output,
        rules.install_rlib_strandphaser.output
    message: 'Performing environment setup...'



def make_log_useful(log_path, status):

    error_buffer = []
    record = 0
    with open(log_path, 'r') as logfile:
        for line in logfile:
            if 'error' in line.lower() or 'exception' in line.lower():
                if record == 0:
                    error_buffer.append('=== Recording ERROR entry')
                record = 10
                error_buffer.append(line.strip())
            elif line.lower().startswith('[w::') or line.lower().startswith('[e::'):
                if record == 0:
                    error_buffer.append('=== Recording library stderr')
                record = 3
                error_buffer.append(line.strip())
            elif record > 0:
                error_buffer.append(line.strip())
                record -= 1
                if record == 0:
                    error_buffer.append('=== STOP recording error entry')
            else:
                continue

    with open(log_path, 'w') as logfile:
        _ = logfile.write('\n'.join(error_buffer))
        _ = logfile.write('\n\n')

    my_env = dict(os.environ)
    with open(log_path, 'a') as logfile:
        _ = logfile.write('\n===[{}]===\n'.format(status))
        _ = logfile.write('Host: {}\n'.format(my_env.get('HOST', 'N/A')))
        _ = logfile.write('Hostname: {}\n'.format(my_env.get('HOSTNAME', 'N/A')))
        _ = logfile.write('Display: {}\n'.format(my_env.get('DISPLAY', 'N/A')))
        _ = logfile.write('Shell: {}\n'.format(my_env.get('SHELL', 'N/A')))
        _ = logfile.write('Terminal: {}\n'.format(my_env.get('TERM', 'N/A')))
        _ = logfile.write('Screen: {}\n'.format(my_env.get('STY', 'N/A')))
        _ = logfile.write('Conda ENV: {}\n'.format(my_env.get('CONDA_DEFAULT_ENV', 'N/A')))
        _ = logfile.write('\n')
    return


onsuccess:
    make_log_useful(log, 'SUCCESS')
    if config['notify']:
        shell('mail -s "[Snakemake] DGA - SUCCESS" {} < {{log}}'.format(config['notify_email']))


onerror:
    make_log_useful(log, 'ERROR')
    if config['notify']:
        shell('mail -s "[Snakemake] DGA - ERRROR" {} < {{log}}'.format(config['notify_email']))
