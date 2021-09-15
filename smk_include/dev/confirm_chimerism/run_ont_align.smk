
localrules: dump_contig_name, run_all

rule run_all:
    input:
        'output/ont_align/NA18989_ONTUL_MAP-TO_chimeric-contigs.mmap.paf.gz'


rule dump_contig_name:
    """
    contig names from file
    2021_pgas/run_folder/output/reference_assembly/clustered/temp/saarclust/results/NA18989_hgsvc_pbsq2-ccs_1000_nhr-hifiasm/NA18989_hgsvc_ilnxs-80pe_sseq/clustered_assembly/
    asmErrorsReport_2e+05bp_dynamic.tsv
    ~5 Mbp
    ~8 Mbp
    ~13 Mbp
    """
    output:
        'output/chimeric.contigs.txt'
    run:
        contigs = [
            'ptg000187l',
            'ptg000086l',
            'ptg000068l'
        ]
        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(contigs) + '\n')


rule extract_contig_sequences:
    input:
        names = 'output/chimeric.contigs.txt',
        assembly = '/gpfs/project/projects/medbioinf/projects/hgsvc/2021_pgas/run_folder/output/reference_assembly/non-hap-res/NA18989_hgsvc_pbsq2-ccs_1000_nhr-hifiasm.fasta'
    output:
        'output/NA1898.chimeric-contigs.fasta'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    shell:
        'seqtk subseq {input.assembly} {input.names} > {output}'


rule align_ont_reads:
    """
    """
    input:
        reads = '/gpfs/project/projects/medbioinf/projects/hgsvc/ontqc/run_folder/input/ONTUL/NA18989_ONTUL_guppy-5.0.11-sup-prom.fasta.gz',
        reference = 'output/NA1898.chimeric-contigs.fasta',
    output:
        paf = 'output/ont_align/NA18989_ONTUL_MAP-TO_chimeric-contigs.mmap.paf.gz',
    benchmark:
        'rsrc/output/ont_align/NA18989_ONTUL_MAP-TO_chimeric-contigs.mmap.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 24
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: 2 * attempt,
    params:
        preset = lambda wildcards: set_alignment_preset(wildcards),
    shell:
        'minimap2 -t 22 -x map-ont -k17 --secondary=no '
        '--cap-kalloc=1g -K4g '
        '{input.reference} {input.reads} | pigz --best -p 4 > {output}'