
include: 'targets.smk'

localrules: run_pur_pr05_family,
          run_pur_pr05_child,
          run_pur_pr05_father,
          run_pur_pr05_mother


rule run_pur_pr05_child:
    input:
         'output/targets/PUR_PR05/HG00733.fofn'
    message: 'Creating results for PUR family PR05 child (HG00733)'


rule run_pur_pr05_father:
    input:
         'output/targets/PUR_PR05/HG00731.fofn'
    message: 'Creating results for PUR family PR05 father (HG00731)'


rule run_pur_pr05_mother:
    input:
         'output/targets/PUR_PR05/HG00732.fofn'
    message: 'Creating results for PUR family PR05 mother (HG00732)'


rule run_pur_pr05_family:
    input:
         rules.run_pur_pr05_father.input,
         rules.run_pur_pr05_mother.input,
         rules.run_pur_pr05_child.input
    message: 'Creating results for PUR family PR05 (admixed population)'