
localrules: run_amr_trios,
            run_pur_trio,
            run_pur_father,
            run_pur_mother,
            run_pur_child,
            run_clm_trio,
            run_clm_child,
            run_mxl_trio,
            run_mxl_child,
            run_pel_trio,
            run_pel_child


rule run_pur_father:
    input:
        'output/targets/AMR_PUR_PR05/HG00731.fofn'
    message: 'Running AMR-PUR-PR05 father'


rule run_pur_mother:
    input:
        'output/targets/AMR_PUR_PR05/HG00732.fofn'
    message: 'Running AMR-PUR-PR05 mother'


rule run_pur_child:
    input:
        'output/targets/AMR_PUR_PR05/HG00733.fofn'
    message: 'Running AMR-PUR-PR05 child'


rule run_pur_trio:
    input:
        rules.run_pur_father.input,
        rules.run_pur_mother.input,
        rules.run_pur_child.input
    message: 'Running AMR-PUR-PR05 trio'

#############################################

rule run_clm_child:
    input:
         'output/targets/AMR_CLM_CLM03/HG01114.fofn'
    message: 'Running AMR-CLM-CLM03 child'

rule run_clm_trio:
    input:
         rules.run_clm_child.input
    message: 'Running AMR-CLM-CLM03 trio'

##############################################

rule run_mxl_child:
    input:
         'output/targets/AMR_MXL_m001/NA19650.fofn'
    message: 'Running AMR-MXL-m001 child'

rule run_mxl_trio:
    input:
         rules.run_mxl_child.input
    message: 'Running AMR-MXL-m001 trio'

##############################################

rule run_pel_child:
    input:
         'output/targets/AMR_PEL_PEL003/HG01573.fofn'
    message: 'Running AMR-PEL-PEL003 child'

rule run_pel_trio:
    input:
         rules.run_pel_child.input
    message: 'Running AMR-PEL-PEL003 trio'

##############################################

rule run_amr_trios:
    input:
        rules.run_pur_trio.input,
        rules.run_clm_trio.input,
        rules.run_mxl_trio.input,
        rules.run_pel_trio.input
    message: 'Running AMR trios'

