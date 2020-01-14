
localrules: run_amr_trios,
            run_pur_trio,
            run_pur_father,
            run_pur_mother,
            run_pur_child


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


rule run_amr_trios:
    input:
        rules.run_pur_trio.input
    message: 'Running AMR trios'
