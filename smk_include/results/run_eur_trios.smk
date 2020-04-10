
localrules: run_hg00171_individual,
            run_na20509_individual,
            run_hg00096_individual,
            run_eur_trios,
            run_ceu_trio,
            run_ceu_child,
            run_ibs_trio,
            run_ibs_child,


rule run_hg00171_individual:
    input:
        'output/targets/EUR_FIN_HG00171/HG00171.fofn'
    message: 'Running EUR-FIN-HG00171 individual'

#######################################################

rule run_na20509_individual:
    input:
        'output/targets/EUR_TSI_NA20509/NA20509.fofn'
    message: 'Running EUR-TSI-NA20509 individual'

#######################################################

rule run_hg00096_individual:
    input:
        'output/targets/EUR_GBR_HG00096/HG00096.fofn'
    message: 'Running EUR-GBR-HG00096 individual'

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

rule run_ibs_child:
    input:
        'output/targets/EUR_IBS_IBS002/HG01505.fofn'
    message: 'Running EUR-IBS-IBS002 child'

rule run_ibs_trio:
    input:
        rules.run_ibs_child.input,
    message: 'Running EUR-IBS-IBS002 trio'

########################################################

rule run_eur_trios:
    input:
        rules.run_hg00171_individual.input,
        rules.run_na20509_individual.input,
        rules.run_hg00096_individual.input,
        rules.run_ceu_trio.input,
        rules.run_ibs_trio.input,
    message: 'Running EUR trios'
