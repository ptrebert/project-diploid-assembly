
include: 'prep_custom_references.smk'
include: 'run_kmer_analysis.smk'
include: 'run_illumina_qv.smk'
include: 'run_tech_comparison.smk'
include: 'run_contig_remap.smk'
include: 'run_bng_hybrids.smk'


localrules: master_eval

wildcard_constraints:
    folder_path = '[A-Za-z0-9\-_\/]+',  # note: "." is NOT allowed in a folder path
    file_name = '[A-Za-z0-9\-_\.]+',
    known_ref = 'GRCh3[78][A-Za-z0-9_]+',
    genemodel = 'GRCh38[A-Za-z0-9_]+'


def quast_busco_determine_targets(wildcards):
    """
    Rerun QUAST-LG with manually fixed BUSCO database
    to get BUSCO stats as requested
    ODB source:
    https://busco-data.ezlab.org/v4/data/lineages/eukaryota_odb10.2020-09-10.tar.gz

    NB: this requires a manual fix for QUAST/BUSCO, i.e. adding the above database
    to the correct path in QUAST
    """

    genemodel = 'GRCh38_GENCODEv31_basic'
    refgenome = 'GRCh38_HGSVC2_noalt'
    folder_path = 'evaluation/phased_assemblies'

    fixed_wildcards = {
        'known_ref': refgenome,
        'genemodel': genemodel,
        'folder_path': folder_path
    }

    # output/{folder_path}/{file_name}.fasta'
    target_path = 'output/evaluation/quastlg_busco/{known_ref}-{genemodel}/{folder_path}/{{file_name}}/report.pdf'.format(**fixed_wildcards)

    load_path = os.path.join('output', folder_path)

    phased_assemblies = sorted([f for f in os.listdir(load_path) if f.endswith('.fasta')])

    compute_targets = []
    for ps_assm in phased_assemblies:
        target_file = target_path.format(**{'file_name': ps_assm.strip('.fasta')})
        compute_targets.append(target_file)

    return compute_targets


rule master_eval:
    input:
        tech_comparison_determine_targets,
        kmer_analysis_determine_targets,
        illumina_qv_determine_targets,
        contig_remap_determine_targets,
        bng_hybrids_determine_targets,
        quast_busco_determine_targets