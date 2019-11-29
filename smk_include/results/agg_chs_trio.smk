
include: '../result_vars.smk'
include: 'chs_father_hg00512.smk'
include: 'chs_mother_hg00513.smk'
include: 'chs_child_hg00514.smk'

localrules: run_chs_trio

rule run_chs_trio:
    input:
        rules.run_chs_father.input,
        rules.run_chs_mother.input,
        rules.run_chs_child.input
    message: 'Creating results for CHS family trio (low diversity population)'