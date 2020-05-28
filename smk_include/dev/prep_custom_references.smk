
include: '../constraints.smk'
include: '../aux_utilities.smk'
include: '../handle_reference_download.smk'

rule intersect_highconf_files:
    input:
        hg001 = 'references/downloads/HG001_hg38_giab_highconf.bed',
        hg002 = 'references/downloads/HG002_hg38_giab_highconf.bed',
        hg005 = 'references/downloads/HG005_hg38_giab_highconf.bed',
    output:
        'references/hg38_giab_highconf.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
         'bedtools multiinter -names HG001 HG002 HG005 -i {input} | egrep "HG001,HG002,HG005" > {output}'