
localrules: run_hg01596_individual,
            run_na18534_individual,
            run_na18939_individual,
            run_hg00864_individual,
            run_eas_trios,
            run_chs_trio,
            run_chs_father,
            run_chs_mother,
            run_chs_child,
            run_khv_trio,
            run_khv_child


rule run_hg01596_individual:
    input:
        'output/targets/EAS_KHV_HG01596/HG01596.fofn'
    message: 'Running EAS-KHV-HG01596 individual'

#####################################################

rule run_na18534_individual:
    input:
        'output/targets/EAS_CHB_NA18534/NA18534.fofn'
    message: 'Running EAS-CHB-NA18534 individual'

#####################################################

rule run_na18939_individual:
    input:
        'output/targets/EAS_JPT_NA18939/NA18939.fofn'
    message: 'Running EAS-JPT-NA18939 individual'

#####################################################

rule run_hg00864_individual:
    input:
        'output/targets/EAS_CDX_HG00864/HG00864.fofn'
    message: 'Running EAS-CDX-HG00864 individual'

#####################################################

rule run_chs_mother:
    input:
        'output/targets/EAS_CHS_SH032/HG00512.fofn'
    message: 'Running EAS-CHS-SH032 mother'


rule run_chs_father:
    input:
        'output/targets/EAS_CHS_SH032/HG00513.fofn'
    message: 'Running EAS-CHS-SH032 father'


rule run_chs_child:
    input:
        'output/targets/EAS_CHS_SH032/HG00514.fofn'
    message: 'Running EAS-CHS-SH032 child'


rule run_chs_trio:
    input:
        rules.run_chs_father.input,
        rules.run_chs_mother.input,
        rules.run_chs_child.input
    message: 'Running EAS-CHS-SH032 trio'

#######################################################

rule run_khv_child:
    input:
        'output/targets/EAS_KHV_VN047/HG02018.fofn'
    message: 'Running EAS-KHV-VN047 child'


rule run_khv_trio:
    input:
        rules.run_khv_child.input
    message: 'Running EAS-KHV-VN047 trio'

########################################################

rule run_eas_trios:
    input:
        rules.run_hg01596_individual.input,
        rules.run_na18534_individual.input,
        rules.run_na18939_individual.input,
        rules.run_hg00864_individual.input,
        rules.run_chs_trio.input,
        rules.run_khv_trio.input
    message: 'Running EAS trios'
