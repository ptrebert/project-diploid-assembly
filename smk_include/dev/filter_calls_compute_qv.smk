
import os
import math
import itertools

workdir: '/scratch/bioinf/projects/diploid-genome-assembly/pebert/assembly_qv/run_folder'


URL_RECORDS = {
    'HG001': 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed',
    'HG005': 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/latest/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.bed',
    'HG002': 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'
}


rule master_qv_estimate:
    input:
        'references/hg38_giab_highconf.bed',

        # trio call sets
        expand('output/variant_calls/20-typefilter/{child}_{parent1}_{parent2}_{reads}_aln-to_{child}_{assembly}_hap{haps}.{var_type}.vcf.stats',
                child='HG00733',
                parent1='HG00731',
                parent2='HG00732',
                reads=['short'],
                assembly=['ont', 'clr', 'ccs', 'zev'],
                haps=[1, 2],
                var_type=['snv', 'indels']),
        expand('output/variant_calls/30-gtfilter/{child}_{parent1}_{parent2}_{reads}_aln-to_{child}_{assembly}_hap{haps}.{var_type}.{genotype}.vcf.stats',
                child='HG00733',
                parent1='HG00731',
                parent2='HG00732',
                reads=['short'],
                assembly=['ont', 'clr', 'ccs', 'zev'],
                haps=[1, 2],
                var_type=['snv', 'indels'],
                genotype=['hom', 'het']),
        expand('output/qv_estimates/{individual}_{parent1}_{parent2}_{reads}_aln-to_{individual}_{assembly}_hap{haps}.qv_stats',
               individual='HG00733',
               parent1='HG00731',
               parent2='HG00732',
               reads=['short'],
               assembly=['zev', 'clr', 'ont', 'ccs'],
               haps=[1, 2]),

        # individual call sets
        expand('output/variant_calls/20-typefilter/{individual}_{reads}_aln-to_{individual}_{assembly}_hap{haps}.{var_type}.vcf.stats',
                individual='NA12878',
                reads=['short'],
                assembly=['ccs', 'whd'],
                haps=[1, 2],
                var_type=['snv', 'indels']),
        expand('output/variant_calls/30-gtfilter/{individual}_{reads}_aln-to_{individual}_{assembly}_hap{haps}.{var_type}.{genotype}.vcf.stats',
               individual='NA12878',
               reads=['short'],
               assembly=['ccs', 'whd'],
               haps=[1, 2],
               var_type=['snv', 'indels'],
               genotype=['hom', 'het']),
        expand('output/variant_calls/60-complete/{individual}_{reads}_aln-to_{individual}_{assembly}_hap{haps}.snv.hg38rev.other.bed',
               individual='NA12878',
               reads=['short'],
               assembly=['ccs', 'whd'],
               haps=[1, 2]),
        expand('output/qv_estimates/{individual}_{reads}_aln-to_{individual}_{assembly}_hap{haps}.qv_stats',
               individual='NA12878',
               reads=['short'],
               assembly=['ccs', 'whd'],
               haps=[1, 2]),



collect_callsets = sorted(filter(lambda x: x.endswith('vcf.bgz'), os.listdir('output/variant_calls/00-raw')))
collect_assemblies = sorted(filter(lambda x: '_hap' in x and x.endswith('.fasta'), os.listdir('references')))
collect_cov_stats = sorted(filter(lambda x: x.endswith('.mdup.sort.cov_stats'), os.listdir('output/alignments')))
print('Available call sets')
for c in collect_callsets:
    print(c)

rule produce_precomputed_input:
    input:
        callsets = ancient(sorted([os.path.join('output/variant_calls/00-raw', x) for x in collect_callsets])),
        assemblies = ancient(sorted([os.path.join('references', x) for x in collect_assemblies])),
        cov_stats = ancient(sorted([os.path.join('output/alignments', x) for x in collect_cov_stats]))
    output:
        sorted([os.path.join('output/variant_calls/linked/calls', x) for x in collect_callsets]),
        sorted([os.path.join('output/variant_calls/linked/assemblies', x) for x in collect_assemblies]),
        sorted([os.path.join('output/variant_calls/linked/stats', x) for x in collect_cov_stats]),
    run:
        for i, o in zip(input, output):
            if os.path.islink(os.path.abspath(o)):
                continue
            os.symlink(os.path.abspath(i), os.path.abspath(o))


rule download_highconf_regions:
    output:
        'references/{sample}_hg38_giab_highconf.bed'
    params:
        url = lambda wildcards: URL_RECORDS[wildcards.sample]
    shell:
        'wget --no-verbose -O {output} {params.url}'


rule intersect_highconf_files:
    input:
        hg001 = 'references/HG001_hg38_giab_highconf.bed',
        hg002 = 'references/HG002_hg38_giab_highconf.bed',
        hg005 = 'references/HG005_hg38_giab_highconf.bed',
    output:
        'references/hg38_giab_highconf.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    priority: 100
    shell:
         'bedtools multiinter -names HG001 HG002 HG005 -i {input} | egrep "HG001,HG002,HG005" > {output}'


rule create_bgzip_tbi_index:
    input:
         ancient('{filepath}.vcf.bgz')
    output:
          '{filepath}.vcf.bgz.tbi'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'bcftools index --tbi {input}'


def load_gw_avg_cov_limit(cov_stats, std_mult):

    with open(cov_stats, 'r') as table:
        _ = table.readline()
        genome_row_cols = table.readline().split()
        gw_mean = float(genome_row_cols[2])
        gw_stddev = float(genome_row_cols[3])
    return int(round(gw_mean + std_mult * gw_stddev))


rule quality_filter_callsets:
    input:
        vcf = ancient('output/variant_calls/00-raw/{reads}_aln-to_{assembly}.vcf.bgz'),
        tbi = 'output/variant_calls/00-raw/{reads}_aln-to_{assembly}.vcf.bgz.tbi',
        vcf_stats = 'output/variant_calls/00-raw/{reads}_aln-to_{assembly}.vcf.stats',
        ref = ancient('references/{assembly}.fasta'),
        cov = ancient('output/alignments/{reads}_aln-to_{assembly}.mdup.sort.cov_stats')
    output:
        'output/variant_calls/10-qfilter/{reads}_aln-to_{assembly}.vcf.bgz'
    conda:
        '../../environment/conda/conda_biotools.yml'
    params:
        cov_limit = lambda wildcards, input: load_gw_avg_cov_limit(input.cov, 3)
    shell:
        'bcftools filter --output-type u --include "QUAL>=10 && INFO/DP<{params.cov_limit}" {input.vcf} | '
        'bcftools norm --fasta-ref {input.ref} --output /dev/stdout --output-type v | '
        'bgzip -c /dev/stdin > {output}'


rule split_callsets_by_type:
    input:
         'output/variant_calls/10-qfilter/{callset}.vcf.bgz',
         'output/variant_calls/10-qfilter/{callset}.vcf.bgz.tbi',
         'output/variant_calls/10-qfilter/{callset}.vcf.stats'
    output:
          'output/variant_calls/20-typefilter/{callset}.snv.vcf.bgz',
          'output/variant_calls/20-typefilter/{callset}.indels.vcf.bgz',
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'bcftools view --types snps --max-alleles 4 '
         '--output-type v --output-file /dev/stdout {input[0]} | '
         'bgzip -c /dev/stdin > {output[0]}'
         ' && '
         'bcftools view --types indels --max-alleles 4 '
         '--output-type v --output-file /dev/stdout {input[0]} | '
         'bgzip -c /dev/stdin > {output[1]}'


rule split_callsets_by_genotype:
    input:
         'output/variant_calls/20-typefilter/{callset}.{var_type}.vcf.bgz',
         'output/variant_calls/20-typefilter/{callset}.{var_type}.vcf.bgz.tbi'
    output:
          'output/variant_calls/30-gtfilter/{callset}.{var_type}.hom.vcf.bgz',
          'output/variant_calls/30-gtfilter/{callset}.{var_type}.het.vcf.bgz',
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'bcftools view --genotype hom '
         '--output-type v --output-file /dev/stdout {input[0]} | '
         'bgzip -c /dev/stdin > {output[0]}'
         ' && '
         'bcftools view --genotype het '
         '--output-type v --output-file /dev/stdout {input[0]} | '
         'bgzip -c /dev/stdin > {output[1]}'


rule compute_callset_stats:
    input:
        ancient('{filepath}.vcf.bgz'),
        '{filepath}.vcf.bgz.tbi'
    output:
        '{filepath}.vcf.stats'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'bcftools stats {input[0]} > {output}'


# bcftools view -m 2 -M 2 --types snps HG00733_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p1/20-qfilter/HG0073F_hgsvc_pbsq2-ccs.Q10.vcf | bcftools query -f '[%SAMPLE=%GT\t]\n' | sort | uniq -c | sort -k1n


rule create_contig_to_reference_alignment:
    """
    Usage note: most paftools do not support input generated with
    the "--eqx" option; the below call is thus different from the
    pipeline version producing SAM/BAM output
    """
    input:
         ref = 'references/hg38.fasta',
         assm = ancient('references/{assembly}.fasta')
    output:
          'output/liftover_paf/{assembly}_lift-to_hg38.paf'
    conda:
         '../../environment/conda/conda_biotools.yml'
    threads: 16
    shell:
         'minimap2 -t {threads} -cx asm20 --cs '
         '--secondary=no -Y -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 '
         '{input.ref} {input.assm} > {output}'


rule dump_vcf_snv_to_bed:
    input:
        'output/variant_calls/20-typefilter/{callset}.snv.vcf.bgz'
    output:
        'output/variant_calls/20-typefilter/{callset}.snv.vcf.bed'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'bgzip -d -c < {input} | vcf2bed --do-not-split --snvs > {output}'


rule dump_vcf_insertions_to_bed:
    input:
         'output/variant_calls/20-typefilter/{callset}.indels.vcf.bgz'
    output:
          'output/variant_calls/20-typefilter/{callset}.ins.vcf.bed'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'bgzip -d -c < {input} | vcf2bed --do-not-split --insertions > {output}'


rule dump_vcf_deletions_to_bed:
    input:
         'output/variant_calls/20-typefilter/{callset}.indels.vcf.bgz'
    output:
          'output/variant_calls/20-typefilter/{callset}.dels.vcf.bed'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'bgzip -d -c < {input} | vcf2bed --do-not-split --deletions > {output}'


rule lift_call_sets_to_reference:
    input:
        bed = 'output/variant_calls/20-typefilter/{reads}_aln-to_{assembly}.{var_type}.vcf.bed',
        paf = 'output/liftover_paf/{assembly}_lift-to_hg38.paf'
    output:
        'output/variant_calls/40-lifted/{reads}_aln-to_{assembly}.{var_type}.hg38.vcf.bed'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'paftools.js liftover -l 10000 {input.paf} {input.bed} > {output}'


rule restrict_calls_to_high_conf_regions:
    input:
        calls = 'output/variant_calls/40-lifted/{reads}_aln-to_{assembly}.{var_type}.hg38.vcf.bed',
        regions = 'references/hg38_giab_highconf.bed'
    output:
        'output/variant_calls/50-highconf/{reads}_aln-to_{assembly}.{var_type}.hg38.hc-in.bed',
        'output/variant_calls/50-highconf/{reads}_aln-to_{assembly}.{var_type}.hg38.hc-out.bed'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools intersect -u -a {input.calls} -b {input.regions} > {output[0]}'
        ' && '
        'bedtools intersect -v -a {input.calls} -b {input.regions} > {output[1]}'


rule reverse_reference_coordinates:
    """
    paftools does not lift additional BED columns
    (= full VCF info); just reverse coordinates
    to simplify intersection with full VCF info
    """
    input:
        'output/variant_calls/50-highconf/{reads}_aln-to_{assembly}.{var_type}.hg38.{location}.bed'
    output:
        'output/variant_calls/50-highconf/{reads}_aln-to_{assembly}.{var_type}.hg38rev.{location}.bed'
    run:
        out_buffer = []
        with open(input[0], 'r') as bed:
            ln = 0
            for line in bed:
                ln += 1
                cols = line.strip().split()
                source_contig, source_start, source_end = cols[3].rsplit('_', 2)
                try:
                    _ = int(source_start)
                    _ = int(source_end)
                except ValueError:
                    # sometimes weird garbage (?) suffix
                    try:
                        source_contig, source_start, source_end = cols[3].rsplit('_', 3)[:-1]
                        _ = int(source_start)
                        _ = int(source_end)
                    except ValueError:
                        raise ValueError('{} / {} / {}'.format(input[0], ln, line.strip()))
                reference_coordinates = '_'.join(cols[:3])
                out_buffer.append((source_contig, source_start, source_end, reference_coordinates))

        out_buffer = sorted(out_buffer, key=lambda x: (x[0], int(x[1]), int(x[2])))
        out_buffer = ['\t'.join(t) for t in out_buffer]

        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(out_buffer) + '\n')


rule intersect_full_vcf_info:
    input:
        original = 'output/variant_calls/20-typefilter/{reads}_aln-to_{assembly}.{var_type}.vcf.bed',
        hc_reg = 'output/variant_calls/50-highconf/{reads}_aln-to_{assembly}.{var_type}.hg38rev.hc-in.bed',
        other_reg = 'output/variant_calls/50-highconf/{reads}_aln-to_{assembly}.{var_type}.hg38rev.hc-out.bed'
    output:
        hc_var = 'output/variant_calls/60-complete/{reads}_aln-to_{assembly}.{var_type}.hg38rev.hc-in.bed',
        non_hc_var = 'output/variant_calls/60-complete/{reads}_aln-to_{assembly}.{var_type}.hg38rev.hc-out.bed',
        other_var = 'output/variant_calls/60-complete/{reads}_aln-to_{assembly}.{var_type}.hg38rev.other.bed'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools intersect -u -a {input.original} -b {input.hc_reg} > {output.hc_var}'
        ' && '
        'bedtools intersect -u -a {input.original} -b {input.other_reg} > {output.non_hc_var}'
        ' && '
        'bedtools intersect -v -a {input.original} -b {input.hc_reg} {input.other_reg} > {output.other_var}'


rule extract_original_sample_order:
    input:
        ancient('output/variant_calls/00-raw/{reads}_aln-to_{assembly}.vcf.bgz'),
    output:
        'output/variant_calls/00-raw/{reads}_aln-to_{assembly}.samples'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'bcftools view -h {input} | egrep "\#CHROM" > {output}'


def comp_qv(num_error_bp, bp_ref=3200000000):
    p = (num_error_bp / (bp_ref * 2))
    try:
        q = -10 * math.log10(p)
    except ValueError:
        return 100
    return int(round(q, 0))


def parse_bed_vcf_row(row, samples):

    assert row['format'].startswith('GT'), 'No genotype field: {} / {}'.format(file_indicator, row)
    genotypes = dict()
    for s in samples:
        gt = row[s].split(':')[0]
        if gt == '.':
            gt = 'N/A'
        genotypes[s] = gt.replace('/', '')

    infos = dict(item.split('=') for item in row['info'].split(';'))
    try:
        infos['LEN'] = int(infos['LEN'])
    except ValueError:
        infos['LEN'] = -1
    return genotypes, infos['TYPE'], infos['LEN']


def determine_trio_error(child, parent1, parent2):

    if child == 'NA':
        if parent1 != 'NA' and parent2 != 'NA':
            call_class = 'AMBIG'
        else:
            call_class = 'SKIP'  # child plus parent missing GT
    else:
        if parent1 != 'NA' and parent2 != 'NA':
            possible_genotypes = set()
            for allele1, allele2 in itertools.product(list(parent1), list(parent2)):
                possible_genotypes.add(allele1 + allele2)
                possible_genotypes.add(allele2 + allele1)
            if child in possible_genotypes:
                call_class = 'OK'
            else:
                call_class = 'ERROR'
        else:
            # this means one parent is GT missing
            call_class = 'OK'
    return call_class


rule compute_qv_estimate:
    input:
        vcf_records = expand('output/variant_calls/60-complete/{{reads}}_aln-to_{{assembly}}.{var_type}.hg38rev.{region_type}.bed',
                             var_type=['snv', 'ins', 'dels'],
                             region_type=['hc-in', 'hc-out', 'other']),
        assembly_index = 'references/{assembly}.fasta.fai',  # not needed anymore...
        sample_order = 'output/variant_calls/00-raw/{reads}_aln-to_{assembly}.samples'
    output:
        'output/qv_estimates/{reads}_aln-to_{assembly}.qv_stats'
    run:
        import collections as col
        import csv
        import operator as op

        descendent_samples = ['NA12878', 'HG00733', 'HG00514', 'NA19240']

        region_types = {
            'hc-in': 'high-confidence',
            'hc-out': 'lifted-not-hc',
            'other': 'not-lifted'
        }

        info_counter = col.Counter()

        with open(input.sample_order, 'r') as foo:
            samples = foo.readline().strip().split()[9:]
        assert len(samples) == 1 or len(samples) == 3, 'Can only handle single or trio: {}'.format(samples)

        if len(samples) == 3:
            descendant = [s for s in samples if s in descendent_samples]
            ancestors = sorted([s for s in samples if s not in descendent_samples])
            # for consistent ordering in output
            trio = descendant + ancestors
            info_counter[('SAMPLES', '-'.join(trio))] = 3
            get_genotypes = op.itemgetter(*tuple(trio))
        else:
            info_counter[('SAMPLES', samples[0])] = 1
            get_genotypes = op.itemgetter(*tuple(samples))

        fieldnames = ['contig', 'start', 'end', 'identifier', 'qual', 'ref', 'alt',
                      'filter', 'info', 'format']
        fieldnames.extend(samples)

        observed_genotypes = set()

        for bedfile in input.vcf_records:
            file_indicator = os.path.basename(bedfile).rsplit('.', 2)[1]
            region_type = region_types[file_indicator]
            with open(bedfile, 'r', newline='') as table:
                ln = 0
                reader = csv.DictReader(table, fieldnames=fieldnames, delimiter='\t')
                for row in reader:
                    ln += 1
                    info_counter[('record', region_type)] += 1
                    genotypes, var_type, var_length = parse_bed_vcf_row(row, samples)
                    [observed_genotypes.add(gt) for gt in genotypes.values()]
                    if var_length < 0:
                        if len(samples) == 1:
                            info_counter[('SKIP', get_genotypes(genotypes), var_type, region_type, 'count')] += 1
                        else:
                            info_counter[('SKIP', '-'.join(get_genotypes(genotypes)), var_type, region_type, 'count')] += 1
                        continue
                    if len(samples) == 1:
                        gt = get_genotypes(genotypes)
                        call_class = 'ERROR' if get_genotypes(genotypes) in ['00', '11'] else 'OK'
                    else:
                        f, p1, p2 = get_genotypes(genotypes)
                        gt = '-'.join([f, p1, p2])
                        if f in ['00', '11']:
                            call_class = 'ERROR'
                        else:
                            call_class = determine_trio_error(f, p1, p2)
                    info_counter[(call_class, gt, var_type, region_type, 'num-bp')] += var_length
                    info_counter[(call_class, gt, var_type, 'any', 'num-bp')] += var_length
                    info_counter[(call_class, gt, 'any', region_type, 'num-bp')] += var_length
                    info_counter[(call_class, gt, 'any', 'any', 'num-bp')] += var_length

        key_sets = {
            ('QV', 'snp', 'error', 'high-conf'): [('ERROR', 'snp', 'high-confidence')],
            ('QV', 'snp', 'error_ambig', 'high-conf'): [('ERROR', 'snp', 'high-confidence'),
                                                        ('AMBIG', 'snp', 'high-confidence')],
            ('QV', 'snp', 'error', 'hc-or-not-lifted'): [('ERROR', 'snp', 'high-confidence'),
                                                         ('ERROR', 'snp', 'not-lifted')],
            ('QV', 'snp', 'error_ambig', 'hc-or-not-lifted'): [('ERROR', 'snp', 'high-confidence'),
                                                               ('AMBIG', 'snp', 'high-confidence'),
                                                               ('ERROR', 'snp', 'not-lifted'),
                                                               ('AMBIG', 'snp', 'not-lifted')],

            ('QV', 'all', 'error', 'high-conf'): [('ERROR', 'high-confidence')],
            ('QV', 'all', 'error_ambig', 'high-conf'): [('ERROR', 'high-confidence'),
                                                        ('AMBIG', 'high-confidence')],
            ('QV', 'all', 'error', 'hc-or-not-lifted'): [('ERROR', 'high-confidence'),
                                                         ('ERROR', 'not-lifted')],
            ('QV', 'all', 'error_ambig', 'hc-or-not-lifted'): [('ERROR', 'high-confidence'),
                                                               ('AMBIG', 'high-confidence'),
                                                               ('ERROR', 'not-lifted'),
                                                               ('AMBIG', 'not-lifted')],
        }

        counter_keys = list(info_counter.keys())
        for label, key_set in key_sets.items():
            if len(samples) == 1 and '_ambig' in label[2]:
                # no ambig calls for single callsets
                continue
            all_selected_keys = []
            for keys in key_set:
                selected_keys = filter(lambda x: all(k in x for k in keys), counter_keys)
                all_selected_keys.extend(selected_keys)
            total_num_bp = sum([info_counter[k] for k in all_selected_keys])
            info_counter[label] = comp_qv(total_num_bp)
            input_num_bp = *label, 'num-bp'
            info_counter[input_num_bp] = total_num_bp

        with open(output[0], 'w') as stats:
            for key in sorted(info_counter.keys()):
                value = info_counter[key]
                _ = stats.write('{}\t{}\n'.format('..'.join(key), value))
            _ = stats.write('observed_genotypes\t{}\n'.format(','.join(sorted(observed_genotypes))))
