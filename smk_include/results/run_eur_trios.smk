
localrules: run_eur_trios,
            run_hg00171_individual,

rule run_hg00171_individual:
    input:
        'output/targets/EUR_FIN_HG00171/HG00171.fofn'
    message: 'Running SAS-ITU-HG00171 individual'


rule run_eur_trios:
    input:
        rules.run_hg00171_individual.input
    message: 'Running EUR trios'
