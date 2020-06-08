
import functools

include: '../constraints.smk'
include: '../aux_utilities.smk'
include: '../handle_reference_download.smk'
include: '../preprocess_references.smk'

wildcard_constraints:
    known_ref = 'GRCh3[78]_[A-Za-z0-9]+_[A-Za-z0-9]+',


rule remove_mitochondrial_sequence:
    input:
        ref_assm = 'references/assemblies/{known_ref}.fasta',
        ref_sizes = 'references/assemblies/{known_ref}.sizes',
    output:
        ref_assm = 'references/assemblies/{known_ref}.no-mito.fasta',
        ref_sizes = 'references/assemblies/{known_ref}.no-mito.sizes',
    resources:
        mem_total_mb = lambda wildcards, attempt: 6144 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 6144 * attempt
    run:
        import io

        buffer = io.StringIO()
        with open(input.ref_assm, 'r') as fasta:
            for line in fasta:
                if line.startswith('>chrM'):
                    skip = True
                elif line.startswith('>'):
                    skip = False
                if skip:
                    continue
                buffer.write(line)

        with open(output.ref_assm, 'w') as dump:
            _ = dump.write(buffer.getvalue())

        buffer = io.StringIO()
        with open(input.ref_sizes, 'r') as sizes:
            for line in sizes:
                if line.startswith('chrM'):
                    continue
                buffer.write(line)

        with open(output.ref_sizes, 'w') as dump:
            dump.write(buffer.getvalue())


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
        'run/references/annotation/{known_ref}-{annotation}.fasta.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt
    shell:
         'bedtools getfasta -fo {output} -name+ -fi {input.fasta} -bed {input.bed}'