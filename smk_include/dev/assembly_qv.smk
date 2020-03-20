
workdir: '/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/run_folder'

WORKDIR = '/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/run_folder'

SHORT_READS_METADATA = '/home/pebert/work/code/github/project-diploid-assembly/annotation/PRJEB9396.tsv'

URL_RECORDS = dict()

rule master_assembly_qv:
    input:
        expand('input/fastq/{individual}_1kg_il25k-125pe_short_{mate}.fastq.gz',
               individual=['HG00733', 'HG00732', 'HG00731'],
               mate=[1, 2]),
        expand('output/alignments/{individual}_short_aln-to_HG00733_hap{num}.mdup.sort.bam{ext}',
               individual=['HG00733', 'HG00732', 'HG00731'],
               num=[1, 2],
               ext=['', '.bai']),
        expand('output/alignments/{individual}_ccs_aln-to_HG00733_hap{num}.mdup.sort.bam{ext}',
               individual=['HG00733', 'HG00732', 'HG00731'],
               num=[1, 2],
               ext=['', '.bai']),


rule link_supp_data:
    input:
        ancient('/scratch/bioinf/users/pebert/data_sources/pur_trio_hifi/HG00731_hgsvc_pbsq2-ccs.fastq.gz'),
        ancient('/scratch/bioinf/users/pebert/data_sources/pur_trio_hifi/HG00731_hgsvc_pbsq2-ccs.fastq.gz'),
        ancient('/scratch/bioinf/users/pebert/data_sources/pur_trio_hifi/HG00731_hgsvc_pbsq2-ccs.fastq.gz'),
        ancient('/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/references/assemblies/hg38_GCA_p13.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200309_Sequel2-CCS_HG00733_V8/HG00733_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200309_Sequel2-CCS_HG00733_V8/HG00733_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200313_Sequel2-CCS_HG00731_HG00732_polished/HG00731_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200313_Sequel2-CCS_HG00731_HG00732_polished/HG00731_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200313_Sequel2-CCS_HG00731_HG00732_polished/HG00732_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta'),
        ancient('/MMCI/TM/scratch/pebert/share/globus/out_hgsvc/20200313_Sequel2-CCS_HG00731_HG00732_polished/HG00732_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta')
    output:
        'input/fastq/HG00731_ccs.fastq.gz',
        'input/fastq/HG00732_ccs.fastq.gz',
        'input/fastq/HG00733_ccs.fastq.gz',
        'references/hg38.fasta',
        'references/HG00733_hap1.fasta',
        'references/HG00733_hap2.fasta',
        'references/HG00731_hap1.fasta',
        'references/HG00731_hap2.fasta',
        'references/HG00732_hap1.fasta',
        'references/HG00732_hap2.fasta',
    run:
        for infile, outfile in zip(input, output):
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
    shell:
        'gzip -d -c {input} | gzip > {output}'


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
         mate1 = 'input/fastq/{individual}_1kg_il25k-125pe_short_1.fastq.gz',
         mate2 = 'input/fastq/{individual}_1kg_il25k-125pe_short_2.fastq.gz',
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
        preset = lambda wildcards: 'CCS' if wildcards.tech == 'ccs' else 'SUBREAD'
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
    threads: 8
    shell:
        'sambamba markdup -t {threads} --overflow-list-size 600000 {input} {output}'


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


# rule align_long_reads_to_haploid_assembly: