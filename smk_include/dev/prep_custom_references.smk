
import functools

include: '../module_includes.smk'

localrules: reduce_to_4col_bed

wildcard_constraints:
    known_ref = 'GRCh3[78]_[A-Za-z0-9]+_[A-Za-z0-9]+',


rule separate_mitochondrial_sequence:
    input:
        ref_assm = 'references/assemblies/{known_ref}.fasta',
        ref_sizes = 'references/assemblies/{known_ref}.sizes',
    output:
        ref_assm = 'references/assemblies/{known_ref}.no-mito.fasta',
        ref_sizes = 'references/assemblies/{known_ref}.no-mito.sizes',
        seq_mito = 'references/assemblies/{known_ref}.chrM.fasta',
    resources:
        mem_total_mb = lambda wildcards, attempt: 6144 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 6144 * attempt
    run:
        import io

        buffer = io.StringIO()
        mito_buffer = io.StringIO()
        with open(input.ref_assm, 'r') as fasta:
            for line in fasta:
                if line.startswith('>chrM'):
                    mito_buffer.write('>chrM\n')
                    skip = True
                    continue
                elif line.startswith('>'):
                    skip = False
                if skip:
                    mito_buffer.write(line)
                    continue
                buffer.write(line)

        with open(output.ref_assm, 'w') as dump:
            _ = dump.write(buffer.getvalue())
        
        with open(output.seq_mito, 'w') as dump:
            _ = dump.write(mito_buffer.getvalue())

        buffer = io.StringIO()
        with open(input.ref_sizes, 'r') as sizes:
            for line in sizes:
                if line.startswith('chrM'):
                    continue
                buffer.write(line)

        with open(output.ref_sizes, 'w') as dump:
            dump.write(buffer.getvalue())


rule retain_alt_seq_in_reference:
    input:
        hgsvc_chroms = 'references/assemblies/GRCh38_HGSVC2_noalt.sizes',
        full_ref = 'references/assemblies/GRCh38_GCA_p13.fasta'
    output:
        ref_with_alt = 'references/assemblies/GRCh38_HGSVC2_incalt.fasta',
        ref_with_alt_sizes = 'references/assemblies/GRCh38_HGSVC2_incalt.sizes'
    resources:
        mem_total_mb = lambda wildcards, attempt: 6144 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 6144 * attempt
    run:
        import io
        import collections

        hgsvc_chroms = set()
        with open(input.hgsvc_chroms, 'r') as table:
            for line in table:
                hgsvc_chroms.add(line.split()[0])

        chrom_sizes = collections.Counter()

        out_buffer = io.StringIO()
        keep_chrom = None
        with open(input.full_ref, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    chrom_name = line.strip('>').strip()
                    if chrom_name.endswith('_alt'):
                        keep_chrom = chrom_name
                    elif chrom_name.endswith('_random'):
                        keep_chrom = chrom_name
                    elif chrom_name.startswith('chrUn'):
                        keep_chrom = chrom_name
                    else:
                        keep_chrom = chrom_name if chrom_name in hgsvc_chroms else None
                    if keep_chrom is not None:
                        out_buffer.write(line)
                    continue
                if keep_chrom is not None:
                    chrom_sizes[chrom_name] += len(line.strip())
                    out_buffer.write(line)
        
        with open(output.ref_with_alt, 'w') as dump:
            _ = dump.write(out_buffer.getvalue())

        with open(output.ref_with_alt_sizes, 'w') as dump:
            for name, length in chrom_sizes.most_common():
                _ = dump.write('{}\t{}\n'.format(name, length))
    # END OF RUN BLOCK


rule intersect_highconf_files:
    input:
        hg001 = 'references/downloads/HG001_hg38_giab_highconf.bed',
        hg002 = 'references/downloads/HG002_hg38_giab_highconf.bed',
        hg005 = 'references/downloads/HG005_hg38_giab_highconf.bed',
    output:
        'references/annotation/hg38_giab_highconf.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
         'bedtools multiinter -names HG001 HG002 HG005 -i {input} | egrep "HG001,HG002,HG005" > {output}'


@functools.lru_cache(maxsize=4096)
def match_chromosome_names(references, query):

    import difflib

    sm = difflib.SequenceMatcher(autojunk=True)
    sm.set_seq2(query)
    len_query = len(query)
    longest_match = 0
    matched_chrom = None
    for refseq in references:
        sm.set_seq1(refseq)
        match = sm.find_longest_match(0, len(refseq), 0, len_query)
        if match.size > longest_match and match.size > 3:
            longest_match = match.size
            matched_chrom = refseq
    return matched_chrom


rule gff_to_bed_match_chromosome_names:
    input:
        gff = 'references/downloads/{ref_assm}_{annotation}.gff3.gz',
        ref_chroms = 'references/assemblies/{ref_assm}_{assm_version}.no-mito.sizes'
    output:
        bedfile = 'references/annotation/{ref_assm}_{assm_version}-{annotation}.bed',
        statsfile = 'references/annotation/{ref_assm}_{assm_version}-{annotation}.stats'
    params:
        min_region_length = 31
    wildcard_constraints:
        ref_assm = 'GRCh[0-9]+',
        assm_version = '[A-Z0-9a-z]+_[A-Z0-9a-z]+'
    run:
        import gzip
        import csv
        import re
        import collections as col

        extract_info = 'ENS[A-Z0-9]+'

        header = ['chrom', 'annotation', 'region_type', 'start', 'end',
                  'score', 'strand', 'frame', 'attribute']

        type_maps = {
            'CTCF_binding_site': 'CTCF',
            'enhancer': 'ENH',
            'open_chromatin_region': 'OPEN',
            'promoter': 'PROM',
            'promoter_flanking_region': 'PRFL',
            'TF_binding_site': 'TFBS'
        }

        with open(input.ref_chroms, 'r') as table:
            ref_chrom_names = set([l.split()[0] for l in table.readlines()])

        match_chromosomes = functools.partial(match_chromosome_names, tuple(ref_chrom_names))

        uniq_match = col.defaultdict(set)

        buffers_rows = []
        kept_entries = col.Counter()
        kept_bp = col.Counter()
        skipped_entries = col.Counter()
        skipped_bp = col.Counter()
        with gzip.open(input.gff, 'rt', newline='') as annotation:
            reader = csv.DictReader(annotation, delimiter='\t', fieldnames=header)
            for row in reader:
                regtype = type_maps.get(row['region_type'], row['region_type'])
                s = int(row['start']) - 1  # GFF is 1-based
                e = int(row['end'])
                if e - s < params.min_region_length:
                    skipped_entries[regtype] += 1
                    skipped_bp[regtype] += e - s
                    continue
                row['start'] = s
                row['end'] = e
                if row['chrom'] in ref_chrom_names:
                    pass
                elif 'chr' + row['chrom'] in ref_chrom_names:
                    row['chrom'] = 'chr' + row['chrom']
                else:
                    matched_ref = match_chromosomes(row['chrom'])
                    if matched_ref is None:
                        skipped_entries[row['chrom']] += 1
                        skipped_entries[regtype] += 1
                        skipped_bp[row['chrom']] += e - s
                        skipped_bp[regtype] += e - s
                        continue
                    uniq_match[matched_ref].add(row['chrom'])
                    row['chrom'] = matched_ref
                kept_entries[regtype] += 1
                kept_bp[regtype] += e - s
                kept_entries['TOTAL'] += 1
                kept_bp['TOTAL'] += e - s
                kept_entries[row['chrom']] += 1
                kept_bp[row['chrom']] += e - s

                add_info = ''
                mobj = re.search(extract_info, row['attribute'])
                if mobj is not None:
                    add_info = mobj.group()
                if add_info and len(add_info) > 3:
                    regtype += '_{}'.format(add_info)
                row['name'] = regtype
                buffers_rows.append(row)

        buffers_rows = sorted(buffers_rows, key=lambda d: (d['chrom'], d['start'], d['end']))

        non_uniq_matches = [(k, v) for k, v in uniq_match.items() if len(v) > 1]
        if non_uniq_matches:
            raise ValueError('Non-unique chromosome matching: {}'.format(non_uniq_matches))

        with open(output.bedfile, 'w', newline='') as bedfile:
            _ = bedfile.write('#')
            writer = csv.DictWriter(bedfile,
                                    fieldnames=['chrom', 'start', 'end', 'name', 'score', 'strand'],
                                    delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(buffers_rows)

        with open(output.statsfile, 'w') as table:
            _ = table.write('# RETAINED\n')
            for key in sorted(kept_entries.keys()):
                _ = table.write('{}_count\t{}\n'.format(key, kept_entries[key]))
                _ = table.write('{}_bp\t{}\n'.format(key, kept_bp[key]))
            _ = table.write('\n# DISCARDED\n')
            for key in sorted(skipped_entries.keys()):
                _ = table.write('{}_count\t{}\n'.format(key, skipped_entries[key]))
                _ = table.write('{}_bp\t{}\n'.format(key, skipped_bp[key]))

            _ = table.write('\n# CHROM-MATCHING\n')
            for key in sorted(uniq_match.keys()):
                _ = table.write('{}\t{}\n'.format(key, uniq_match[key].pop()))
    ### END OF RULE


rule add_sequences_to_bed:
    input:
        bed = 'references/annotation/{known_ref}-{annotation}.bed',
        fasta = 'references/assemblies/{known_ref}.no-mito.fasta',
        fai = 'references/assemblies/{known_ref}.no-mito.fasta.fai'
    output:
        'references/annotation/{known_ref}-{annotation}.fasta'
    benchmark:
        'rsrc/references/annotation/{known_ref}-{annotation}.fasta.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt
    shell:
         'bedtools getfasta -fo {output} -name+ -fi {input.fasta} -bed {input.bed}'


rule profile_sequences_in_bed:
    """
    Limit BED input, some characters in the input
    may crash bedtools nuc / core dump or exception
    """
    input:
        bed = 'references/annotation/{annotation}.4c.bed',
        fasta = 'references/assemblies/{known_ref}.fasta',
        fai = 'references/assemblies/{known_ref}.fasta.fai'
    output:
        'references/annotation/{known_ref}_{annotation}.nuc.stats'
    log:
        'log/references/annotation/{known_ref}_{annotation}.nuc.stats.log'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt
    shell:
         'bedtools nuc -fi {input.fasta} -bed {input.bed} > {output} 2> {log}'


rule convert_fasta_to_hdf:
    input:
        'references/annotation/{known_ref}-{annotation}.fasta'
    output:
        'references/annotation/{known_ref}-{annotation}.h5'
    benchmark:
        'rsrc/references/annotation/{known_ref}-{annotation}.hdf.rsrc'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd

        seq_df = []
        with open(input[0], 'r') as fasta:
            seq_buffer = ''
            for line in fasta:
                if line.startswith('>'):
                    if seq_buffer:
                        seq_df.append((regtype, regid, chrom, int(start), int(end), seq_buffer))
                        seq_buffer = ''
                    regtype, rest = line[1:].strip().split('_', 1)
                    regid, _, chrom, coords = rest.split(':')
                    start, end = coords.split('-')
                else:
                    seq_buffer += line.strip()
            seq_df.append((regtype, regid, chrom, int(start), int(end), seq_buffer))

        seq_df = pd.DataFrame.from_records(
            seq_df,
            columns=[
                'region_type', 'region_id', 'chrom',
                'start', 'end', 'hg38_seq'
            ])
        seq_df.to_hdf(output[0], 'sequences')
    # end of run block


rule gzip_delta_file:
    input:
        'output/evaluation/{filepath}.delta'
    output:
        'output/evaluation/{filepath}.delta.gz'
    shell:
        'gzip -c {input} > {output}'


rule gunzip_trio_assembly:
    input:
        'references/downloads/{assembly}.fa.gz'
    output:
        'references/assemblies/{assembly}.fasta'
    wildcard_constraints:
        assembly = '[HGNA]{2}[0-9]{5}_hgsvc_pbil-trio_hap[AB]'
    shell:
        'gzip -d -c {input} > {output}'


rule reduce_to_4col_bed:
    """
    bigWigAverageOverBed is somewhat picky and seems to check
    if the last column (not the 6th...) is a valid strand
    identifier ('+' or '-') in case there are more than 4
    columns in the BED file
    """
    input:
        '{filepath}/{filename}.bed'
    output:
        '{filepath}/{filename}.4c.bed'
    shell:
        'cut -f 1,2,3,4 {input} > {output}'
    

rule prep_pav_calls:
    """
    PAV call set here is not the published one (i.e. Zenodo)
    variant IDs may no longer match with published version,
    in particular relevant for the list of dropped variants
    used here
    MD5: 36c95ab5684535033bcea988ceb8c3e3
    """
    input:
        good = 'references/downloads/variants_freeze{version}_sv_insdel.prepublish.tsv.gz',
        bad = 'references/downloads/variants-dropped_freeze{version}_sv_insdel.tsv.gz'
    output:
        bad = 'references/annotation/PAV_sv-insdel-dropped_v{version}.bed',
        both = 'references/annotation/PAV_sv-insdel_v{version}.h5'
    benchmark:
        'rsrc/references/annotation/PAV_sv-insdel_v{version}.hdf.rsrc'
    run:
        import pandas as pd

        lowq = pd.read_csv(input.bad, sep='\t', index_col=False, encoding='ascii', dtype=str)
        lowq['quality'] = '0'

        hiq = pd.read_csv(input.good, sep='\t', index_col=False, encoding='ascii', dtype=str)
        hiq['quality'] = '1'

        bed_out = lowq.loc[:, ['#CHROM', 'POS', 'END', 'ID']].copy()
        bed_out['POS'] = bed_out['POS'].astype('int64')
        bed_out['END'] = bed_out['END'].astype('int64')
        bed_out.sort_values(['#CHROM', 'POS'], inplace=True)
        bed_out.to_csv(output.bad, sep='\t', index=False, header=False)

        both = pd.concat([hiq, lowq], axis=0)
        both['POS'] = both['POS'].astype('int64')
        both['END'] = both['END'].astype('int64')
        both.sort_values(['#CHROM', 'POS', 'END'], inplace=True)
        both['POS'] = both['POS'].astype(str)
        both['END'] = both['END'].astype(str)
        both.reset_index(drop=True, inplace=True)

        both.to_hdf(output.both, key='PAV_v{}'.format(wildcards.version), mode='w', format='fixed')
