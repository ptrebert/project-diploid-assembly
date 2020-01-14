
localrules: run_eas_trios,
            run_chs_trio,
            run_chs_father,
            run_chs_mother,
            run_chs_child


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


rule run_eas_trios:
    input:
        rules.run_chs_trio.input
    message: 'Running EAS trios'
