localrules: master

samples = [
    'HG00512',
    'HG00514',
    'HG00731',
    'HG00732',
    'HG00733',
    'NA19238',
    'NA19239',
    'NA19240'
]

output_files = []

for hap in ['H1', 'H2']:
    for s in samples:
        tmp = 'output/tigs/{sample}_{hap}_tigs.fasta'.format(**{'sample': s, 'hap': hap})
        output_files.append(tmp)

rule extract_hap1_contigs:
    input:
        tig_names = '/beeond/data/hifiasm_v13_assemblies/tig_names/{sample}_H1_tigs.txt',
        assm_fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h1-un.fasta',
    output:
        assm_tigs = 'output/tigs/{sample}_H1_tigs.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'seqtk subseq -l 120 {input.assm_fasta} {input.tig_names} > {output}'


rule extract_hap2_contigs:
    input:
        tig_names = '/beeond/data/hifiasm_v13_assemblies/tig_names/{sample}_H2_tigs.txt',
        assm_fasta = '/beeond/data/hifiasm_v13_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-hifiasm.h2-un.fasta',
    output:
        assm_tigs = 'output/tigs/{sample}_H2_tigs.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'seqtk subseq -l 120 {input.assm_fasta} {input.tig_names} > {output}'
        
rule master:
    input:
        output_files
