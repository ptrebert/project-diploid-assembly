
localrules: master_statistics_input_data

rule master_statistics_input_data:
    input:
        []


rule compute_statistics_complete_input_fastq:
    input:
        fastq = 'input/fastq/{sample}.fastq.gz',
        faidx = 'references/assemblies/' + config['use_genome_size'] +'.fasta.fai',
    output:
        dump = 'output/statistics/stat_dumps/{sample}.fastq.pck',
        summary = 'input/fastq/{sample}.stats',
    log: 'log/output/statistics/stat_dumps/{sample}.fastq.log',
    benchmark: 'rsrc/output/statistics/stat_dumps/{sample}.fastq.t2.rsrc'
    threads: 2
    resources:
        runtime_hrs= 8,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 2048,
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.fastq} '
        '--output {output.dump} --summary-output {output.summary} '
        '--copy-stats-dump '
        'output/statistics/stat_dumps/{wildcards.sample}.pbn.bam.pck '
        'output/statistics/stat_dumps/{wildcards.sample}.fasta.pck '
        '--copy-summary '
        'input/bam/{wildcards.sample}.stats '
        'input/fasta/{wildcards.sample}.stats '
        '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule compute_statistics_split_cluster_fastq:
    input:
         fasta = 'output/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.gz',
         faidx = 'output/reference_assembly/clustered/{sseq_reads}/{reference}/sequences/{sequence}.seq',
    output:
          dump = 'output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.pck',
          summary = 'output/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.stats',
    log: 'log/output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.log',
    benchmark: 'rsrc/output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.t2.rsrc'
    threads: 2
    resources:
             runtime_hrs= 1,
             mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
             mem_per_cpu_mb = lambda wildcards, attempt: (1024 + 1024 * attempt) // 2,
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
          script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
         '{params.script_exec} --debug --input-files {input.fasta} '
         '--output {output.dump} --summary-output {output.summary} '
         '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule compute_statistics_joint_cluster_fastq:
    input:
        fasta = 'output/diploid_assembly/strandseq_joint/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.fastq.gz',
        faidx = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta.fai',
    output:
        dump = 'output/statistics/stat_dumps/diploid_assembly/strandseq_joint/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.fastq.pck',
        summary = 'output/diploid_assembly/strandseq_joint/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.stats',
    log:
        'log/output/statistics/stat_dumps/diploid_assembly/strandseq_joint/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.fastq.log',
    benchmark:
        'rsrc/output/statistics/stat_dumps/diploid_assembly/strandseq_joint/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/draft/haploid_fastq/{hap_reads}.{hap}.fastq.t2.rsrc'
    threads: 2
    resources:
        runtime_hrs= lambda wildcards, attempt: attempt,
        mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: (1024 + 1024 * attempt) // 2,
    conda:
        '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.fasta} '
        '--output {output.dump} --summary-output {output.summary} '
        '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule compute_statistics_complete_input_bam:
    input:
        bam = 'input/bam/{sample}.pbn.bam',
        faidx = 'references/assemblies/' + config['use_genome_size'] +'.fasta.fai',
    output:
        dump = 'output/statistics/stat_dumps/{sample}.pbn.bam.pck',
        summary = 'input/bam/{sample}.stats',
    log: 'log/output/statistics/stat_dumps/{sample}.pbn.bam.log',
    benchmark: 'rsrc/output/statistics/stat_dumps/{sample}.pbn.bam.t2.rsrc'
    threads: 2
    resources:
        runtime_hrs= 23,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 2048,
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.bam} '
        '--output {output.dump} --summary-output {output.summary} '
        '--copy-stats-dump '
        'output/statistics/stat_dumps/{wildcards.sample}.fasta.pck '
        'output/statistics/stat_dumps/{wildcards.sample}.fastq.pck '
        '--copy-summary '
        'input/fasta/{wildcards.sample}.stats '
        'input/fastq/{wildcards.sample}.stats '
        '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule collect_snv_stats_per_cluster:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf',
        fai = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta.fai'
    output:
        'output/statistics/variant_calls/{var_caller}/{reference}/{sseq_reads}/{vc_reads}.snv.QUAL{qual}.GQ{gq}.vcf.cluster.stats'
    priority: 200
    run:
        import statistics as stats
        import collections as col

        with open(input.fai, 'r') as index:
            cluster_sizes = dict((l.split()[0], int(l.split()[1])) for l in index.readlines())

        snv_per_chrom = col.Counter()
        qual_per_snv = col.defaultdict(list)
        depth_per_snv = col.defaultdict(list)
        genoqual_per_snv = col.defaultdict(list)
        with open(input.vcf, 'r') as vcf:
            for line in vcf:
                if line.startswith('#'):
                    continue
                chrom, _, _, _, _, qual, _, info, attributes, sample = line.strip().split()
                snv_per_chrom[chrom] += 1
                qual_per_snv[chrom].append(float(qual))
                if wildcards.var_caller == 'longshot':
                    assert info.startswith('DP='), 'Unexpected Longshot/INFO field (DP): {}'.format(line.strip())
                    depth_per_snv[chrom].append(int(info.split(';')[0].replace('DP=', '')))
                elif wildcards.var_caller == 'deepvar' or wildcards.var_caller == 'freebayes':
                    assert attributes.split(':')[2] == 'DP', 'Unexpected DeepVariant/FORMAT field (DP): {}'.format(line.strip())
                    depth_per_snv[chrom].append(int(sample.split(':')[2]))
                else:
                    raise RuntimeError('Unsupported variant caller for cluster stats rule: {}'.format(wildcards.var_caller))
                assert attributes.split(':')[1] == 'GQ', 'Unexpected FORMAT field (GQ): {}'.format(line.strip())
                genoqual_per_snv[chrom].append(int(sample.split(':')[1]))

        fun_labels = ['mean', 'stddev']
        stats_funs = [stats.mean, stats.stdev]

        stats_labels = ['QUAL', 'DEPTH', 'GTQUAL']
        stats_collect = [qual_per_snv, depth_per_snv, genoqual_per_snv]

        with open(output[0], 'w') as stat_dump:
            for seq in sorted(cluster_sizes.keys()):
                _ = stat_dump.write('{}_size_bp\t{}\n'.format(seq, cluster_sizes[seq]))
                _ = stat_dump.write('{}_HET-SNV_num\t{}\n'.format(seq, snv_per_chrom[seq]))
                kbp_factor = cluster_sizes[seq] / 1000
                snv_per_kbp = round(snv_per_chrom[seq] / kbp_factor, 3)
                _ = stat_dump.write('{}_HET-SNV_per_kbp\t{}\n'.format(seq, snv_per_kbp))
                for sl, sc in zip(stats_labels, stats_collect):
                    for fl, fun in zip(fun_labels, stats_funs):
                        try:
                            result = round(fun(sc[seq]), 3)
                        except stats.StatisticsError:
                            # can happen for no data
                            result = 0.
                        _ = stat_dump.write('{}_HET-SNV_{}_{}\t{}\n'.format(seq, sl, fl, result))


def select_assembly_contig_index(wildcards):

    formatter = dict(wildcards)
    regular_fasta_idx = 'output/{folder_path}/{assembly}.fasta.fai'.format(**formatter)
    split_fasta_idx = regular_fasta_idx.replace('.fasta.fai', '.fasta-split.fa.fai')
    if 'clustered' in split_fasta_idx:
        return split_fasta_idx
    else:
        return regular_fasta_idx


rule collect_contig_to_ref_aln_statistics:
    input:
        ctg_ref_aln = 'output/alignments/contigs_to_reference/{folder_path}/{assembly}_map-to_{aln_reference}.bed',
        ref_chroms = 'references/assemblies/{aln_reference}.fasta.fai',
        assm_chroms = select_assembly_contig_index
    output:
        'output/statistics/contigs_to_ref_aln/{folder_path}/{assembly}_map-to_{aln_reference}.mapq{mapq}.stats'
    priority: 200
    conda:
        '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('collect_contig_aln_stats.py'),
        mapq = lambda wildcards: int(wildcards.mapq),
        chrom_cov = lambda wildcards: int(config['min_chromosome_coverage'][wildcards.aln_reference]),
        no_fail = '--ignore-min-chrom-cov' if config['ignore_min_chrom_cov'] else ''
    shell:
        '{params.script_exec} --contig-alignments {input.ctg_ref_aln} '
        '--reference-chromosomes {input.ref_chroms} --contig-names {input.assm_chroms} '
        '--min-chrom-coverage {params.chrom_cov} {params.no_fail} '
        '--min-mapq {params.mapq} --contig-groups --group-id-position 0 --output {output}'


rule compute_assembly_contig_summary:
    input:
        assm = 'output/{folder_path}/{assembly}.fasta.fai',
        ref = 'references/assemblies/' + config['use_genome_size'] +'.fasta.fai',
    output:
        'output/statistics/assembly_summary/{folder_path}/{assembly}.stats.tsv'
    priority: 200
    conda:
        '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('collect_contig_stats.py'),
    shell:
        '{params.script_exec} --assembly-fai {input.assm} --reference-fai {input.ref} --output {output} '


def collect_tag_lists(wildcards, glob_collect=True, caller='snakemake'):
    """
    :param wildcards:
    :return:
    """
    import os

    source_path = os.path.join('output',
                               PATH_STRANDSEQ_DGA_SPLIT,
                               'draft/haplotags/{hap_reads}.{sequence}.tags.{tag_type}.tsv')
    glob_path = source_path.replace('{sequence}', '*')
    cluster_key = 'sequence'
    if glob_collect:
        seq_files = check_cluster_file_completeness(wildcards, source_path, glob_path, wildcards.sseq_reads, cluster_key)
        if not seq_files:
            raise RuntimeError('collect_tag_lists: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used!')

        reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sseq_reads)
        seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                        reference=wildcards.reference).output[0]
        checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.seq'))

        seq_files = expand(source_path,
                           var_caller=wildcards.var_caller,
                           gq=wildcards.gq,
                           qual=wildcards.qual,
                           reference=wildcards.reference,
                           vc_reads=wildcards.vc_reads,
                           sseq_reads=wildcards.sseq_reads,
                           hap_reads=wildcards.hap_reads,
                           sequence=checkpoint_wildcards.sequence,
                           tag_type=wildcards.tag_type)

    assert seq_files, 'collect_tag_lists >> returned empty output'
    return seq_files


rule summarize_tagging_splitting_statistics:
    """
    2021-05-22 see rule
    integrative_phasing::write_strandphaser_split_vcf_fofn
    for details about fofn input
    """
    input:
        fofn = 'output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.fofn',
        tags = collect_tag_lists
    output:
        'output/statistics/tag_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{hap_reads}.tags.{tag_type}.tsv'
    benchmark:
        'rsrc/output/statistics/tag_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{hap_reads}.tags.{tag_type}.rsrc'
    priority: 200
    run:
        import os
        import collections as col

        # Sanity check: there must be one tag list per cluster
        sample_name = wildcards.sseq_reads.split('_')[0]
        num_clusters = estimate_number_of_saarclusters(sample_name, wildcards.sseq_reads, readset=wildcards.hap_reads)

        num_tags = len(input.tags)

        num_wh_vcf = 0
        with open(input.fofn, 'r') as fofn:
            num_wh_vcf = len([line for line in fofn if line.strip()])
        assert num_wh_vcf > 0, 'write_assembled_fasta_clusters_fofn >> number of WhatsHap VCF splits read from fofn is zero: {}'.format(input.fofn)

        if num_tags == 0:
            raise RuntimeError('summarize_tagging_splitting_statistics >> zero haplo-tag files: {}'.format(wildcards))
        elif num_tags != num_clusters and num_tags != num_wh_vcf:
            raise RuntimeError('summarize_tagging_splitting_statistics >> mismatch between expected ({}) and received ({}) haplotag files: {}'.format(num_clusters, num_tags, wildcards))
        else:
            pass

        hapcount = col.Counter()
        line_num = 0

        for tsv in sorted(input.tags):
            with open(tsv, 'r') as table:
                for line in table:
                    if line.startswith('#'):
                        continue
                    line_num += 1
                    try:
                        hap = line.split()[1]
                    except IndexError:
                        raise IndexError('Malformed line {} from file {}: "{}"'.format(line_num, os.path.basename(tsv), line))
                    hapcount[hap] += 1

        entry_types = tuple(sorted(hapcount.keys()))
        if entry_types != ('H1', 'H2', 'none'):
            raise ValueError('Haplotype entries in tag list not as expected (H1/H2/none): {}'.format(entry_types))

        total_reads = sum(list(hapcount.values()))

        percent_h1 = round((hapcount['H1'] / total_reads) * 100, 2)
        percent_h2 = round((hapcount['H2'] / total_reads) * 100, 2)
        percent_tag = round(percent_h1 + percent_h2, 2)
        percent_un = round((hapcount['none'] / total_reads) * 100, 2)

        with open(output[0], 'w') as summary:
            _ = summary.write('num_reads\t{}\n'.format(total_reads))
            _ = summary.write('num_hap1\t{}\n'.format(hapcount['H1']))
            _ = summary.write('num_hap2\t{}\n'.format(hapcount['H2']))
            _ = summary.write('num_untagged\t{}\n'.format(hapcount['none']))
            _ = summary.write('percent_hap1\t{}\n'.format(percent_h1))
            _ = summary.write('percent_hap2\t{}\n'.format(percent_h2))
            _ = summary.write('percent_tagged\t{}\n'.format(percent_tag))
            _ = summary.write('percent_untagged\t{}\n'.format(percent_un))
        
