
localrules: run_afr_trios,
            run_yri_trio,
            run_yri_father,
            run_yri_mother,
            run_yri_child,
            run_gwd_trio,
            run_gwd_child


rule run_yri_mother:
    input:
        'output/targets/AMR_YRI_Y117/NA19238.fofn'
    message: 'Running AFR-YRI-Y117 mother'


rule run_yri_father:
    input:
        'output/targets/AMR_YRI_Y117/NA19239.fofn'
    message: 'Running AFR-YRI-Y117 father'


rule run_yri_child:
    input:
        'output/targets/AMR_YRI_Y117/NA19240.fofn'
    message: 'Running AFR-YRI-Y117 child'


rule run_yri_trio:
    input:
        rules.run_yri_father.input,
        rules.run_yri_mother.input,
        rules.run_yri_child.input
    message: 'Running AFR-YRI-Y117 trio'

#####################################################

rule run_gwd_child:
    input:
        'output/targets/AFR_GWD_GB24/HG02587.fofn'
    message: 'Running AFR-GWD-GB24 child'


rule run_gwd_trio:
    input:
        rules.run_gwd_child.input,
    message: 'Running AFR-GWD-GB24 trio'

######################################################

rule run_afr_trios:
    input:
        rules.run_yri_trio.input,
        rules.run_gwd_trio.input
    message: 'Running AFR trios'
