
localrules: run_sas_trios,
            run_itu_trio,
            run_itu_mother,

rule run_itu_mother:
    input:
        'output/targets/SAS_ITU_SH032/HG03721.fofn'
    message: 'Running SAS-ITU-IT003 mother'


rule run_itu_trio:
    input:
        rules.run_itu_mother.input,
    message: 'Running SAS-ITU-IT003 trio'


rule run_sas_trios:
    input:
        rules.run_itu_trio.input
    message: 'Running SAS trios'
