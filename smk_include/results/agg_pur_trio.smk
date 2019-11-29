
include: '../result_vars.smk'
include: 'pur_father_hg00731.smk'
include: 'pur_mother_hg00732.smk'
include: 'pur_child_hg00733.smk'

localrules: run_pur_trio

rule run_pur_trio:
    input:
        rules.run_pur_father.input,
        rules.run_pur_mother.input,
        rules.run_pur_child.input
    message: 'Creating results for PUR family trio (admixed population)'