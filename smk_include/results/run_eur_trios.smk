
localrules: run_eur_trios,
            run_hg00171_individual,
            run_ceu_trio,
            run_ceu_child

rule run_hg00171_individual:
    input:
        'output/targets/EUR_FIN_HG00171/HG00171.fofn'
    message: 'Running EUR-FIN-HG00171 individual'

#######################################################

rule run_ceu_child:
    input:
        'output/targets/EUR_CEU_1328/NA12329.fofn'
    message: 'Running EUR-CEU-1328 child'

rule run_ceu_trio:
    input:
        rules.run_ceu_child.input,
    message: 'Running EUR-CEU-1328 trio'

########################################################

rule run_eur_trios:
    input:
        rules.run_hg00171_individual.input,
        rules.run_ceu_trio.input,
    message: 'Running EUR trios'
