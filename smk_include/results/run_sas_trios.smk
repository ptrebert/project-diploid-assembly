
localrules: run_hg03009_individual,
            run_na20847_individual,
            run_sas_trios,
            run_itu_trio,
            run_itu_mother,
            run_itu_child,
            run_stu_trio,
            run_stu_child,
            run_pjl_trio,
            run_pjl_child


rule run_hg03009_individual:
    input:
        'output/targets/SAS_BEB_HG03009/HG03009.fofn'
    message: 'Running SAS-BEB-HG03009 individual'

########################################################

rule run_na20847_individual:
    input:
        'output/targets/SAS_GIH_NA20847/NA20847.fofn'
    message: 'Running SAS-GIH-NA20847 individual'

########################################################

rule run_itu_mother:
    input:
        'output/targets/SAS_ITU_IT003/HG03721.fofn'
    message: 'Running SAS-ITU-IT003 mother'

rule run_itu_child:
    input:
        'output/targets/SAS_ITU_IT003/HG03732.fofn'
    message: 'Running SAS-ITU-IT003 child'

rule run_itu_trio:
    input:
        rules.run_itu_mother.input,
        rules.run_itu_child.input
    message: 'Running SAS-ITU-IT003 trio'

#########################################################

rule run_stu_child:
    input:
        'output/targets/SAS_STU_ST012/HG03683.fofn'
    message: 'Running SAS-STU-ST012 child'

rule run_stu_trio:
    input:
        rules.run_stu_child.input
    message: 'Running SAS-STU-ST012 trio'

#########################################################

rule run_pjl_child:
    input:
        'output/targets/SAS_PJL_PK06/HG02492.fofn'
    message: 'Running SAS-PJL-PK06 child'

rule run_pjl_trio:
    input:
        rules.run_pjl_child.input
    message: 'Running SAS-PJL-PK06 trio'

#########################################################

rule run_sas_trios:
    input:
        rules.run_hg03009_individual.input,
        rules.run_na20847_individual.input,
        rules.run_itu_trio.input,
        rules.run_stu_trio.input,
        rules.run_pjl_trio.input
    message: 'Running SAS trios'
