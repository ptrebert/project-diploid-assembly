
workdir: '/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/run_folder'

WORKDIR = '/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/run_folder'

SHORT_READS_METADATA = '/home/pebert/work/code/github/project-diploid-assembly/annotation/PRJEB9396.tsv'

URL_RECORDS = dict()

rule master_assembly_qv:
    input:
        []
    message: 'Manual target specification is required'

rule link_supp_data:
    input:
        ancient('/scratch/bioinf/users/pebert/data_sources/pur_trio_hifi/HG00731_hgsvc_pbsq2-ccs.fastq.gz'),
        ancient('/scratch/bioinf/users/pebert/data_sources/pur_trio_hifi/HG00731_hgsvc_pbsq2-ccs.fastq.gz'),
        ancient('/scratch/bioinf/users/pebert/data_sources/pur_trio_hifi/HG00731_hgsvc_pbsq2-ccs.fastq.gz'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/references/assemblies/hg38_GCA_p13.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200330_Sequel2-CCS_HG00733_v9_data-collection/haploid_assemblies/HG00733_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200330_Sequel2-CCS_HG00733_v9_data-collection/haploid_assemblies/HG00733_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200313_Sequel2-CCS_HG00731_HG00732_polished/HG00731_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200313_Sequel2-CCS_HG00731_HG00732_polished/HG00731_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200313_Sequel2-CCS_HG00731_HG00732_polished/HG00732_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200313_Sequel2-CCS_HG00731_HG00732_polished/HG00732_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/HG00512_hgsvc_pbsq2-ccs_1000.fastq.gz'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/HG00513_hgsvc_pbsq2-ccs_1000.fastq.gz'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/HG00514_hgsvc_pbsq2-ccs_1000.fastq.gz'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/NA19238_hgsvc_pbsq2-ccs_1000.fastq.gz'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/NA19239_hgsvc_pbsq2-ccs_1000.fastq.gz'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/NA19240_hgsvc_pbsq2-ccs_1000.fastq.gz'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200318_ONT_HG00733_assemblies/HG00733_hpg_ontpm-ul_1000-flye.h1-un.helpol-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200318_ONT_HG00733_assemblies/HG00733_hpg_ontpm-ul_1000-flye.h2-un.helpol-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200309_Sequel2-CLR_HG00512_polished/HG00512_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200309_Sequel2-CLR_HG00512_polished/HG00512_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200311_Sequel2-CLR_HG00513_polished/HG00513_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200311_Sequel2-CLR_HG00513_polished/HG00513_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200302_Sequel2-CLR_HG00514_polished/HG00514_hgsvc_pbsq2-clr_0526-flye.h1-un.arrow-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200302_Sequel2-CLR_HG00514_polished/HG00514_hgsvc_pbsq2-clr_0526-flye.h2-un.arrow-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200306_Sequel2-CCS_NA12878_NA24385_polished/NA12878_giab_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200306_Sequel2-CCS_NA12878_NA24385_polished/NA12878_giab_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/NA12878_giab_pbsq2-ccs_1000.fastq.gz'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200402_Sequel1-CLR_HG00733/HG00733_sra_pbsq1-clr_1000-flye.h1-un.arrow-p1.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200402_Sequel1-CLR_HG00733/HG00733_sra_pbsq1-clr_1000-flye.h2-un.arrow-p1.fasta'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/ext_data/HG00733_zev_hap1.fasta'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/ext_data/HG00733_zev_hap2.fasta'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/ext_data/HG00733_uwwh_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/ext_data/HG00733_uwwh_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta')
    output:
        'input/fastq/HG00731_ccs.fastq.gz',
        'input/fastq/HG00732_ccs.fastq.gz',
        'input/fastq/HG00733_ccs.fastq.gz',
        'references/hg38.fasta',
        'references/HG00733_ccs_hap1.fasta',
        'references/HG00733_ccs_hap2.fasta',
        'references/HG00731_ccs_hap1.fasta',
        'references/HG00731_ccs_hap2.fasta',
        'references/HG00732_ccs_hap1.fasta',
        'references/HG00732_ccs_hap2.fasta',
        'input/fastq/HG00512_ccs.fastq.gz',
        'input/fastq/HG00513_ccs.fastq.gz',
        'input/fastq/HG00514_ccs.fastq.gz',
        'input/fastq/NA19238_ccs.fastq.gz',
        'input/fastq/NA19239_ccs.fastq.gz',
        'input/fastq/NA19240_ccs.fastq.gz',
        'references/HG00733_ont_hap1.fasta',
        'references/HG00733_ont_hap2.fasta',
        'references/HG00512_clr_hap1.fasta',
        'references/HG00512_clr_hap2.fasta',
        'references/HG00513_clr_hap1.fasta',
        'references/HG00513_clr_hap2.fasta',
        'references/HG00514_clr_hap1.fasta',
        'references/HG00514_clr_hap2.fasta',
        'references/NA12878_ccs_hap1.fasta',
        'references/NA12878_ccs_hap2.fasta',
        'input/fastq/NA12878_ccs.fastq.gz',
        'references/HG00733_clr_hap1.fasta',
        'references/HG00733_clr_hap2.fasta',
        'references/HG00733_zev_hap1.fasta',
        'references/HG00733_zev_hap2.fasta',
        'references/HG00733_wh_hap1.fasta',
        'references/HG00733_wh_hap2.fasta'
    run:
        for infile, outfile in zip(input, output):
            if os.path.islink(outfile):
                continue
            os.symlink(infile, outfile)


def select_short_read_fastq_parts(wildcards):
    import csv
    import subprocess

    individual = wildcards.sample.split('_')[0]
    rel_dir = 'input/fastq/{}_1kg_il25k-125pe_short'.format(individual)
    base_name = '{}_1kg_il25k-125pe_'.format(individual)
    fastq_files = []
    with open(SHORT_READS_METADATA, 'r', newline='') as metadata:
        reader = csv.DictReader(metadata, delimiter='\t')
        for row in reader:
            if individual not in row['sample_alias']:
                continue
            run_id = row['run_accession']
            mate1, mate2 = row['fastq_ftp'].split(';')
            if wildcards.mate == '1':
                rel_local_mate1 = os.path.join(rel_dir, base_name + run_id + '_1.fastq.gz')
                URL_RECORDS[rel_local_mate1] = 'ftp://' + mate1
                fastq_files.append(rel_local_mate1)
            else:
                rel_local_mate2 = os.path.join(rel_dir, base_name + run_id + '_2.fastq.gz')
                URL_RECORDS[rel_local_mate2] = 'ftp://' + mate2
                fastq_files.append(rel_local_mate2)
    if not fastq_files:
        raise ValueError('No input fastq files selected')
    return sorted(fastq_files)


rule download_fastq:
    output:
        'input/fastq/{sample}_short/{sample}_{run_id}_{mate}.fastq.gz'
    params:
        dl_url = lambda wildcards, output: URL_RECORDS[output[0]]
    shell:
        'wget --no-verbose -O {output} {params.dl_url}'


rule merge_illumina_fastq_files:
    input:
        select_short_read_fastq_parts
    output:
        'input/fastq/{sample}_short_{mate}.fastq.gz'
    wildcard_constraints:
        sample = '(HG00731|HG00732|HG00733|HG00512|HG00513|HG00514)'
    shell:
        'gzip -d -c {input} | gzip > {output}'


ceph_family_short_reads_url = {
    ('NA12878', '1'): "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz",  # daughter
    ('NA12878', '2'): "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz",
    ('NA12891', '1'): "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194160/ERR194160_1.fastq.gz",  # father
    ('NA12891', '2'): "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194160/ERR194160_2.fastq.gz",
    ('NA12892', '1'): "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194161/ERR194161_1.fastq.gz",
    ('NA12892', '2'): "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194161/ERR194161_2.fastq.gz"
}


rule load_na12878_data:
    """
    bioproject source: PRJEB3381
    alternative: PRJNA200694
    """
    output:
        'input/fastq/{individual}_short_{mate}.fastq.gz',
    wildcard_constraints:
        individual = '(NA12878|NA12891|NA12892)'
    params:
        url = lambda wildcards: ceph_family_short_reads_url[(wildcards.individual, wildcards.mate)]
    shell:
        'wget --no-verbose -O {output} {params.url}'


rule generate_bwa_index:
    input:
         reference = 'references/{assembly}.fasta'
    output:
          'references/{assembly}/bwa_index/{assembly}.amb',
          'references/{assembly}/bwa_index/{assembly}.ann',
          'references/{assembly}/bwa_index/{assembly}.bwt',
          'references/{assembly}/bwa_index/{assembly}.pac',
          'references/{assembly}/bwa_index/{assembly}.sa'
    conda:
         '../../environment/conda/conda_biotools.yml'
    params:
          prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
         'bwa index -p {params.prefix} {input.reference}'


rule bwa_short_reads_to_reference_alignment:
    input:
         mate1 = 'input/fastq/{individual}_short_1.fastq.gz',
         mate2 = 'input/fastq/{individual}_short_2.fastq.gz',
         ref_idx = 'references/{assembly}/bwa_index/{assembly}.bwt'
    output:
          bam = 'output/alignments/{individual}_short_aln-to_{assembly}.sort.bam'
    conda:
         '../../environment/conda/conda_biotools.yml'
    threads: 48
    params:
          idx_prefix = lambda wildcards, input: input.ref_idx.rsplit('.', 1)[0],
          keep_flag = 2,
          discard_flag = 3840
    shell:
         'bwa mem -t {threads}'
         ' -R "@RG\\tID:{wildcards.individual}_short_{wildcards.assembly}\\tPL:Illumina\\tSM:{wildcards.individual}"'
         ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} | '
         ' samtools view -u -F {params.discard_flag} -f {params.keep_flag} /dev/stdin | '
         ' samtools sort -@ 12 -l 6 -m 16G -O BAM -o {output} /dev/stdin'


rule pbmm2_long_reads_to_reference_alignment:
    input:
         reads = 'input/fastq/{individual}_{tech}.fastq.gz',
         reference = 'references/{assembly}.fasta'
    output:
          bam = 'output/alignments/{individual}_{tech}_aln-to_{assembly}.sort.bam'
    conda:
         '../../environment/conda/conda_pbtools.yml'
    threads: 48
    params:
        tempdir = lambda wildcards: os.path.join('temp', 'pbmm2', wildcards.assembly, wildcards.individual, wildcards.tech),
        preset = lambda wildcards: 'CCS --min-length 5000 ' if wildcards.tech == 'ccs' else 'SUBREAD --min-length 5000 '
    shell:
         'TMPDIR={params.tempdir} '
         'pbmm2 align --log-level INFO --sort --sort-memory 16384M --no-bai '
         ' --alignment-threads 40 --sort-threads 8 --preset {params.preset} '
         ' --rg "@RG\\tID:{wildcards.individual}_{wildcards.tech}_{wildcards.assembly}\\tSM:{wildcards.individual}" '
         ' --sample {wildcards.individual} {input.reference} {input.reads} {output.bam}'


rule mark_duplicate_reads:
    input:
        'output/alignments/{alignment}.sort.bam'
    output:
        'output/alignments/{alignment}.mdup.sort.bam'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: 12
    params:
        tempdir = lambda wildcards: '--tmpdir=/local/tmp' if os.path.isdir('/local/tmp') else ''
    shell:
        'sambamba markdup -t {threads} '
        '--sort-buffer-size=8092 '
        '--overflow-list-size 1000000 '  # default 200 000 ; increase to avoid "too many open files" issue
        '--hash-table-size 524288 '  # default 262 144 ; increase to avoid "too many open files" issue
        '{params.tempdir} '
        '{input} {output}'


rule index_bam_alignment:
    input:
         bam = '{filepath}.bam'
    output:
         bai = '{filepath}.bam.bai'
    threads: 12
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         "samtools index -@ {threads} {input.bam}"


rule index_fasta_reference:
    input:
        'references/{reference}.fasta'
    output:
        'references/{reference}.fasta.fai'
    wildcard_constraints:
        reference = '\w+'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'samtools faidx {input[0]}'


rule create_tbi_index:
    input:
         '{filepath}.vcf.bgz'
    output:
          '{filepath}.vcf.bgz.tbi'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'bcftools index --tbi {input}'


checkpoint generate_sequence_files:
    input:
        'references/{assembly}.fasta.fai'
    output:
        directory('references/{assembly}/sequences/')
    wildcard_constraints:
        assembly = '\w+'
    run:
        import os
        os.makedirs(output[0], exist_ok=True)
        with open(input[0], 'r') as index:
            for line in index:
                seq_file = os.path.join(output[0], line.split()[0] + '.seq')
                with open(seq_file, 'w') as dump:
                    _ = dump.write(line)


rule compute_alignments_stats:
    input:
        bam = 'output/alignments/{reads}_aln-to_{assembly}.mdup.sort.bam',
        bai = 'output/alignments/{reads}_aln-to_{assembly}.mdup.sort.bam.bai',
        seq = 'references/{assembly}/sequences/{seq}.seq'
    output:
        'output/alignments/{reads}_aln-to_{assembly}/rmdup_stats/{seq}.stats'
    conda:
         '../../environment/conda/conda_biotools.yml'
    threads: 2
    shell:
        'samtools stats --remove-dups --threads {threads} {input.bam} {wildcards.seq} > {output}'


def select_sequence_stats_files(wildcards):

    stats_files = expand('output/alignments/{individual}_{tech}_aln-to_{assembly}/rmdup_stats/{seq}.stats',
                         individual = wildcards.individuals.split('_'),
                         tech = wildcards.tech,
                         assembly = wildcards.assembly,
                         seq = wildcards.seq)
    return sorted(stats_files)


def weighted_avg_and_std(values, weights):
    """
    https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    import numpy as np
    try:
        average = np.average(values, weights=weights)
    except ZeroDivisionError:
        return 0, 0
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return average, np.sqrt(variance)


rule compute_coverage_statistics:
    input:
        stats_files = select_sequence_stats_files,
        seq_info = 'references/{assembly}/sequences/{seq}.seq'
    output:
        'output/alignments/{individuals}_{tech}_aln-to_{assembly}/cov_stats/{seq}.cov'
    run:
        import numpy as np

        with open(input.seq_info, 'r') as faidx:
            seqlen = int(faidx.readline().split()[1])

        num_inputs = len(list(input.stats_files))

        # coverage range from [ 0 ... 1001 ]
        coverages = np.zeros(1002, dtype=np.int32)

        for stats_file in input.stats_files:
            with open(stats_file, 'r') as stats:
                for line in stats:
                    if not line.startswith('COV'):
                        continue
                    _, cov_range, cov, count = line.strip().split()
                    if '<' in cov_range:
                        cov = 1001
                    else:
                        cov = int(cov)
                    coverages[cov] += int(count)
        # samtools stats excludes bases with coverage 0
        # so the following is the (average) number of bases
        # covered with at least one read per input alignment
        non_zero_cov_bp = int(round(coverages.sum() / num_inputs, 0))

        # compute number of zero coverage bases
        coverages[0] = seqlen * num_inputs - coverages.sum()

        observed_coverage_values = np.array(coverages > 0, dtype=np.bool)
        cov_values = np.array(range(1002), dtype=np.int32)[observed_coverage_values]
        cov_counts = coverages[observed_coverage_values]
        avg, std = weighted_avg_and_std(cov_values, cov_counts)

        with open(output[0], 'w') as dump:
            _ = dump.write('cov_mean\t{}\n'.format(round(avg, 2)))
            _ = dump.write('cov_std\t{}\n'.format(round(std, 2)))
            _ = dump.write('cov_nonzero\t{}\n'.format(non_zero_cov_bp))
            _ = dump.write('seq_length\t{}\n'.format(seqlen))


def select_cov_stats_files(wildcards):
    collect_path = 'output/alignments/{alignment}/cov_stats/{seq}.cov'

    reads, assembly = wildcards.alignment.split('_aln-to_')

    checkpoint_dir = checkpoints.generate_sequence_files.get(assembly=assembly).output[0]

    checkpoint_wildcards = glob_wildcards(os.path.join(checkpoint_dir, '{seq}.seq'))

    stats_files = expand(collect_path,
                         alignment=wildcards.alignment,
                         seq=checkpoint_wildcards.seq)

    return stats_files


rule summarize_average_coverage:
    input:
        select_cov_stats_files
    output:
        'output/alignments/{alignment}.mdup.sort.cov_stats'
    run:
        import numpy

        out_lines = [('sequence', 'length', 'cov_mean', 'cov_std', 'cov_nonzero'), ('genome_wide', 'placeholder')]
        averages = []
        weights = []
        total_length = 0
        total_nonzero = 0

        for sf in sorted(input):
            contig_name = os.path.basename(sf).split('.')[0]
            with open(sf, 'r') as stats:
                _, mean = stats.readline().strip().split()
                _, std = stats.readline().strip().split()
                _, nz_cov = stats.readline().strip().split()
                _, length = stats.readline().strip().split()
                out_lines.append((contig_name, length, mean, std, nz_cov))
                averages.append(float(mean))
                weights.append(int(length))

                total_length += int(length)
                total_nonzero += int(nz_cov)

        averages = numpy.array(averages, dtype=numpy.float32)
        weights = numpy.array(weights, dtype=numpy.int32)

        gw_avg, gw_std = weighted_avg_and_std(averages, weights)
        out_lines[1] = (
            'genome',
            str(total_length),
            str(round(gw_avg, 2)),
            str(round(gw_std, 2)),
            str(total_nonzero)
        )

        with open(output[0], 'w') as table:
            _ = table.write('#')
            _ = table.write('\n'.join(['\t'.join(record) for record in out_lines]))


def compute_coverage_limit(cov_file, num_samples, sd_mult):

    # coverages for Illumina data varies wildly, and freebayes' runtime
    # seems to be very sensitive to coverage. Set some
    # hard limit no matter what the observed coverage actually is
    BOGUS_LIMIT = 400

    with open(cov_file, 'r') as cov_info:
        _, cov_mean = cov_info.readline().split()
        _, cov_stddev = cov_info.readline().split()

    cov_limit = int(round(float(cov_mean) * num_samples + sd_mult * float(cov_stddev), 0))
    cov_limit = max(min(cov_limit, BOGUS_LIMIT), 1)  # max 1: not sure if freebayes would handle 0 correctly
    return str(cov_limit)


rule freebayes_call_variants_trio:
    input:
        child_bam = 'output/alignments/{child}_{reads}_aln-to_{assembly}.mdup.sort.bam',
        child_bai = 'output/alignments/{child}_{reads}_aln-to_{assembly}.mdup.sort.bam.bai',
        parent1_bam = 'output/alignments/{parent1}_{reads}_aln-to_{assembly}.mdup.sort.bam',
        parent1_bai = 'output/alignments/{parent1}_{reads}_aln-to_{assembly}.mdup.sort.bam.bai',
        parent2_bam = 'output/alignments/{parent2}_{reads}_aln-to_{assembly}.mdup.sort.bam',
        parent2_bai = 'output/alignments/{parent2}_{reads}_aln-to_{assembly}.mdup.sort.bam.bai',
        ref = 'references/{assembly}.fasta',
        seq = 'references/{assembly}/sequences/{seq}.seq',
        cov = 'output/alignments/{child}_{parent1}_{parent2}_{reads}_aln-to_{assembly}/cov_stats/{seq}.cov'
    output:
        'output/variant_calls/00-raw/{child}_{parent1}_{parent2}_{reads}_aln-to_{assembly}/splits/{seq}.vcf.bgz'
    log:
       'log/output/variant_calls/00-raw/{child}_{parent1}_{parent2}_{reads}_aln-to_{assembly}/splits/{seq}.freebayes.log'
    params:
        skip_coverage = lambda wildcards, input: compute_coverage_limit(input.cov, 3, 2)
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'freebayes --bam {input.child_bam} --bam {input.parent1_bam} --bam {input.parent2_bam} '
        '--fasta-reference {input.ref} --region {wildcards.seq} --vcf /dev/stdout '
        '--strict-vcf --use-best-n-alleles 6 --skip-coverage {params.skip_coverage} 2> {log} '
        ' | bgzip -c /dev/stdin > {output}'


rule freebayes_call_variants_single:
    input:
         child_bam = 'output/alignments/{individual}_{reads}_aln-to_{assembly}.mdup.sort.bam',
         child_bai = 'output/alignments/{individual}_{reads}_aln-to_{assembly}.mdup.sort.bam.bai',
         ref = 'references/{assembly}.fasta',
         seq = 'references/{assembly}/sequences/{seq}.seq',
         cov = 'output/alignments/{individual}_{reads}_aln-to_{assembly}/cov_stats/{seq}.cov'
    output:
          'output/variant_calls/00-raw/{individual}_{reads}_aln-to_{assembly}/splits/{seq}.vcf.bgz'
    log:
       'log/output/variant_calls/00-raw/{individual}_{reads}_aln-to_{assembly}/splits/{seq}.freebayes.log'
    wildcard_constraints:
        individual = '[A-Z0-9]+'
    params:
        skip_coverage = lambda wildcards, input: compute_coverage_limit(input.cov, 1, 2)
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'freebayes --bam {input.child_bam} '
         '--fasta-reference {input.ref} --region {wildcards.seq} --vcf /dev/stdout '
         '--strict-vcf --use-best-n-alleles 4 --skip-coverage {params.skip_coverage} 2> {log} '
         ' | bgzip -c /dev/stdin > {output}'


def collect_callset_splits(wildcards):
    collect_path = 'output/variant_calls/00-raw/{callset}/splits/{seq}.vcf.bgz{ext}'

    _, assembly = wildcards.callset.split('_aln-to_')

    seq_dir = checkpoints.generate_sequence_files.get(assembly=assembly).output[0]

    checkpoint_wildcards = glob_wildcards(os.path.join(seq_dir, '{seq}.seq'))

    vcf_splits = expand(collect_path,
                        callset=wildcards.callset,
                        seq=checkpoint_wildcards.seq,
                        ext=['', '.tbi'])
    return vcf_splits


rule write_callset_splits_fofn:
    input:
         collect_callset_splits
    output:
          'output/variant_calls/00-raw/{callset}.fofn'
    run:
        import os
        if len(input) < 1:
            raise ValueError('No splits to merge')
        with open(output[0], 'w') as fofn:
            for file_path in sorted(input):
                if file_path.endswith('.tbi'):
                    continue
                _ = fofn.write(file_path + '\n')


rule merge_callset_splits:
    input:
        fofn = 'output/variant_calls/00-raw/{callset}.fofn',
        stats = 'output/alignments/{callset}.mdup.sort.cov_stats'  # just to trigger creating this
    output:
          protected('output/variant_calls/00-raw/{callset}.vcf.bgz')
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'bcftools concat --output /dev/stdout --output-type v --file-list {input.fofn} | bgzip -c /dev/stdin > {output}'
