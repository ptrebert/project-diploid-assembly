
workdir: "/MMCI/TM/scratch/pebert/diploid_assembly/side_tracks/friday_call_set"

chromosomes = ['chr' + str(i) for i in range(1, 23)]
chromosomes.extend(['chrX', 'chrY'])
notify = True


# Email 2019-04-11 Kishwar Shafin shafin@ucsc.edu
friday_call_set = 'https://storage.googleapis.com/kishwar-friday-data/benchmarking/FRIDAY/pacbio_ccs_vcf_friday/pacbio_ccs_hg002_grch38_model_illumina%2Bpacbio.vcf.gz'

# GIAB call sets
call_set_son = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz'
call_set_mother = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz'
call_set_father = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh38/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz'

dl_urls = {'giab_HG002': call_set_son,
            'giab_HG003': call_set_father,
            'giab_HG004': call_set_mother,
            'friday_HG002': friday_call_set}

onsuccess:
    body_text = "Nothing to report\n{log}"
    if notify:
        shell('mail -s "[Snakemake] FRIDAY re-phasing - success" {} <<< "{}"'.format('pebert@mpi-inf.mpg.de', body_text))

onerror:
    if notify:
        shell('mail -s "[Snakemake] FRIDAY re-phasing - ERRROR" {} < {{log}}'.format('pebert@mpi-inf.mpg.de'))


rule main:
    input:
        expand('output/variant_processing/genotyping/GRCh38_friday-hiconf_HG002_{project}.vcf.gz',
                project=['pangen', 'ucsc1', 'pangen-ucsc1'])


rule download_call_sets:
    output:
        'references/GRCh38_{callset}_variants.vcf.gz'
    params:
        dl_url = lambda w: dl_urls[w.callset]
    shell:
        "wget --quiet -O {output} {params.dl_url}"


rule generate_pedigree_file:
    output:
        'references/AJ.ped'
    run:
        with open(output[0], 'w') as pedfile:
            _ = pedfile.write("AJ\tHG002\tHG003\tHG004\t1\t1")


rule build_vcf_tbi_index:
    input:
        '{filepath}.vcf.gz'
    output:
        '{filepath}.vcf.gz.tbi'
    shell:
        'bcftools index --tbi --output-file {output} {input}'


rule enforce_vcf_specification:
    input:
        raw_calls = 'references/GRCh38_friday_HG002_variants.vcf.gz',
        sample_name = 'references/friday_samples.txt'
    output:
        norm = 'references/GRCh38_norm_HG002_variants.vcf.gz',
        reheader = temp('references/friday_reheader.vcf.gz')
    log: 'log/enforce_vcf_spec.log'
    run:
        exec = 'bcftools reheader'
        exec += ' --samples {input.sample_name}'
        exec += ' --output {output.reheader}'
        exec += ' {input.raw_calls}'
        exec += ' && '
        exec += ' bcftools annotate --remove FORMAT'
        exec += ' --output {output.norm} --output-type z'
        exec += ' {output.reheader}'
        exec += ' &> {log}'
        shell(exec)


rule split_friday_calls_by_chrom:
    input:
        call_set = 'references/GRCh38_norm_HG002_variants.vcf.gz',
        call_set_index = 'references/GRCh38_norm_HG002_variants.vcf.gz.tbi'
    output:
        'output/variant_processing/genotyping/preprocessing/GRCh38_norm_HG002_variants.{chrom}.vcf.gz'
    log: 'log/variant_processing/genotyping/preprocessing/GRCh38_norm_HG002_variants.{chrom}.log'
    run:
        exec = 'bcftools view --regions {wildcards.chrom}'
        exec += ' --output-type z --output-file {output}'
        exec += ' {input} &> {log}'
        shell(exec)


rule single_genotype_friday_call_set:
    input:
        callset = 'output/variant_processing/genotyping/preprocessing/GRCh38_norm_HG002_variants.{chrom}.vcf.gz',
        callset_index = 'output/variant_processing/genotyping/preprocessing/GRCh38_norm_HG002_variants.{chrom}.vcf.gz.tbi',
        refgenome = '../../references/GRCh38_decoy_hla.fa',
        refgenome_index = '../../references/GRCh38_decoy_hla.fa.fai',
        alignment = '../../output/alignments/HG002.{project}.ont-ul.sort.cram',
        alignment_index = '../../output/alignments/HG002.{project}.ont-ul.sort.cram.crai',
    output:
        'output/variant_processing/genotyping/GRCh38_norm_HG002_{project}.{chrom}.gt.vcf.gz'
    log: 'log/variant_processing/genotyping/GRCh38_norm_HG002_{project}.{chrom}.gt.log'
    wildcard_constraints:
        project = '[a-z0-9]+',
        chrom = 'chr[0-9XY]+'
    conda:
        "../environment/conda/wh_genotype.yml"
    shell:
        "whatshap --debug genotype --reference {input.refgenome} --output {output} --ignore-read-groups --chromosome {wildcards.chrom} {input.callset} {input.alignment} &> {log}"

#    run:
#        exec = 'whatshap --debug genotype --reference {input.refgenome}'
#        exec += ' --output {output} --ignore-read-groups'
#        exec += ' --chromosome {wildcards.chrom}'
#        exec += ' {input.callset} {input.alignment}'
#        exec += ' &> {log}'
#        shell(exec)


rule multi_genotype_friday_call_set:
    input:
        callset = 'output/variant_processing/genotyping/preprocessing/GRCh38_norm_HG002_variants.{chrom}.vcf.gz',
        callset_index = 'output/variant_processing/genotyping/preprocessing/GRCh38_norm_HG002_variants.{chrom}.vcf.gz.tbi',
        refgenome = '../../references/GRCh38_decoy_hla.fa',
        refgenome_index = '../../references/GRCh38_decoy_hla.fa.fai',
        alignment = expand('../../output/alignments/HG002.{project}.ont-ul.sort.cram',
                            project=['pangen', 'ucsc1']),
        alignment_index = expand('../../output/alignments/HG002.{project}.ont-ul.sort.cram.crai',
                            project=['pangen', 'ucsc1']),
    output:
        'output/variant_processing/genotyping/GRCh38_norm_HG002_pangen-ucsc1.{chrom}.gt.vcf.gz'
    log: 'log/variant_processing/genotyping/GRCh38_norm_HG002_pangen-ucsc1.{chrom}.gt.log'
    wildcard_constraints:
        chrom = 'chr[0-9XY]+'
    conda:
        "../environment/conda/wh_genotype.yml"
    shell:
        "whatshap --debug genotype --reference {input.refgenome} --output {output} --ignore-read-groups --chromosome {wildcards.chrom} {input.callset} {input.alignment} &> {log}"

#    run:
#        exec = 'whatshap --debug genotype --reference {input.refgenome}'
#        exec += ' --output {output} --ignore-read-groups'
#        exec += ' --chromosome {wildcards.chrom}'
#        exec += ' {input.callset} {input.alignment}'
#        exec += ' &> {log}'
#        shell(exec)


rule merge_chromosome_vcf_files:
    input:
        expand('output/variant_processing/genotyping/GRCh38_norm_HG002_{{project}}.{chrom}.gt.vcf.gz',
                chrom=chromosomes)
    output:
        'output/variant_processing/genotyping/GRCh38_norm_HG002_{project}.wg.gt.vcf.gz'
    log: 'log/variant_processing/genotyping/GRCh38_norm_HG002_{project}.merge.log'
    shell:
        'bcftools concat -o {output} -O z {input} &> {log}'


rule intersect_raw_new_genotype_calls:
    input:
        raw = 'references/GRCh38_norm_HG002_variants.vcf.gz',
        raw_index = 'references/GRCh38_norm_HG002_variants.vcf.gz.tbi',
        new = 'output/variant_processing/genotyping/GRCh38_norm_HG002_{project}.wg.gt.vcf.gz',
        new_index = 'output/variant_processing/genotyping/GRCh38_norm_HG002_{project}.wg.gt.vcf.gz.tbi',
    output:
        'output/variant_processing/genotyping/GRCh38_friday-hiconf_HG002_{project}.vcf.gz'
    log: 'log/variant_processing/genotyping/GRCh38_friday-hiconf_HG002_{project}.isect.log'
    run:
        exec = "bcftools isec --include 'TYPE=\"snp\" && GT=\"het\"'"
        exec += " --collapse all -n =2 --output-type z"
        exec += " --output {output}"
        exec += " {input.raw} {input.new}"
        exec += " &> {log}"
        shell(exec)