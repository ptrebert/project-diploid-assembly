
localrules: run_afr_trios,
            run_yri_trio,
            run_yri_father,
            run_yri_mother,
            run_yri_child,
            run_gwd_trio,
            run_gwd_child,
            run_acb_trio,
            run_acb_child,
            run_asw_trio,
            run_asw_child,
            run_msl_trio,
            run_msl_child,
            run_esn_trio,
            run_esn_child

rule run_na19036_individual:
    input:
        'output/targets/AFR_LWK_NA19036/NA19036.fofn'
    message: 'Running AFR-LWK-NA19036 individual'

######################################################

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

rule run_acb_child:
    input:
        'output/targets/AFR_ACB_BB13/HG02011.fofn'
    message: 'Running AFR-ACB-BB13 child'


rule run_acb_trio:
    input:
        rules.run_acb_child.input,
    message: 'Running AFR-ACB-BB13 trio'

######################################################

rule run_asw_child:
    input:
        'output/targets/AFR_ASW_2436/NA19983.fofn'
    message: 'Running AFR-ASW-2436 child'


rule run_asw_trio:
    input:
        rules.run_asw_child.input,
    message: 'Running AFR-ASW-2436 trio'

######################################################

rule run_msl_child:
    input:
        'output/targets/AFR_MSL_SL05/HG03065.fofn'
    message: 'Running AFR-MSL-SL05 child'


rule run_msl_trio:
    input:
        rules.run_asw_child.input,
    message: 'Running AFR-MSL-SL05 trio'

######################################################

rule run_esn_child:
    input:
        'output/targets/AFR_ESN_NG98/HG03371.fofn'
    message: 'Running AFR-ESN-NG98 child'


rule run_esn_trio:
    input:
        rules.run_esn_child.input,
    message: 'Running AFR-ESN-NG98 trio'

######################################################

rule run_afr_trios:
    input:
        rules.run_na19036_individual.input,
        rules.run_yri_trio.input,
        rules.run_gwd_trio.input,
        rules.run_acb_trio.input,
        rules.run_asw_trio.input,
        rules.run_msl_trio.input,
        rules.run_esn_trio.input
    message: 'Running AFR trios'
