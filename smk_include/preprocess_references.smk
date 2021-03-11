
localrules: master_preprocess_references

rule master_preprocess_references:
    input:
        []


rule normalize_reference_assembly_names:
    input:
        'references/downloads/{known_ref}.fa.gz'
    output:
        seq = 'references/assemblies/{known_ref}.fasta',
    resources:
        mem_total_mb = 12288,
        mem_per_cpu_mb = 12288
    run:
        import gzip as gzip
        import re
        import io as io

        get_chrom = re.compile('chromosome\s[0-9XY]+(\s|,)')
        get_patch = re.compile('contig\s[A-Z0-9_\-]+,')

        remove_chars = ['*', ':', '-', '|']

        out_buffer = io.StringIO()
        chrom_length = 0
        current_chrom = None
        with gzip.open(input[0], 'rt') as genome:
            for line in genome:
                if line.startswith('>'):
                    if current_chrom is not None:
                        chrom_length = 0
                    chrom_name = line.split()[0].strip('>')
                    try:
                        int(chrom_name)
                        chrom_name = 'chr' + chrom_name
                    except ValueError:
                        if chrom_name in ['X', 'Y', 'M', 'MT']:
                            chrom_name = 'chr' + chrom_name
                    # some heuristics for nicer and more
                    # consistent naming across all reference
                    # assemblies, in particular for GenBank
                    chrom_name = chrom_name.split('.')[0]
                    if 'mitochondrion' in line.lower() and 'complete' in line.lower():
                        chrom_name = 'chrM'
                    elif 'unplaced genomic contig' in line.lower():
                        chrom_name = 'chrUn_' + chrom_name
                    elif 'unlocalized genomic contig' in line.lower():
                        mobj = get_chrom.search(line)
                        if mobj is not None:
                            chrom_id = mobj.group(0).strip(' ,').split()[1]
                            chrom_name = 'chr' + chrom_id + '_' + chrom_name + '_random'
                    elif 'alternate locus' in line.lower():
                        mobj = get_chrom.search(line)
                        if mobj is not None:
                            chrom_id = mobj.group(0).strip(' ,').split()[1]
                            chrom_name = 'chr' + chrom_id + '_' + chrom_name + '_alt'
                    elif 'PATCH' in line:
                        mobj = get_chrom.search(line)
                        if mobj is not None:
                            chrom_id = mobj.group(0).strip(' ,').split()[1]
                            mobj = get_patch.search(line)
                            if mobj is not None:
                                patch_id = mobj.group(0).strip(' ,').split()[1]
                                patch_id = patch_id.strip('_PATCH')
                                chrom_name = 'chr' + chrom_id + '_' + chrom_name + '_' + patch_id
                    elif 'primary' in line.lower():
                        mobj = get_chrom.search(line)
                        if mobj is not None:
                            chrom_id = mobj.group(0).strip().strip(',').split()[1]
                            chrom_name = 'chr' + chrom_id
                    else:
                        # additional heuristics to be added
                        pass
                    if any([c in chrom_name for c in remove_chars]):
                        for r in remove_chars:
                            chrom_name = chrom_name.replace(r, '_')

                    # end of heuristics
                    current_chrom = chrom_name
                    out_buffer.write('>{}\n'.format(chrom_name))
                else:
                    out_buffer.write(line)
                    chrom_length += len(line.strip())

        with open(output.seq, 'w') as dump:
            _ = dump.write(out_buffer.getvalue())

    # end of rule


rule reduce_reference_to_main_chromosomes:
    input:
        full_ref = 'references/assemblies/GRCh3{num}_{ref_id}.fasta'
    output:
        red_ref = 'references/assemblies/hg3{num}_{ref_id}.fasta',
    wildcard_constraints:
        num = '[0-9]'
    resources:
        mem_total_mb = 12288,
        mem_per_cpu_mb = 12288
    run:
        chrom_set = config['main_chromosomes']

        import io

        out_buffer = io.StringIO()
        skip = True
        current_chrom = None
        chrom_length = 0

        with open(input.full_ref, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    if current_chrom is not None:
                        chrom_length = 0
                    chrom_name = line.strip('>').split()[0]
                    if chrom_name not in chrom_set:
                        skip = True
                        current_chrom = None
                        chrom_length = 0
                        continue
                    skip = False
                    current_chrom = chrom_name
                    out_buffer.write('>{}\n'.format(chrom_name))
                elif skip:
                    continue
                else:
                    out_buffer.write(line)
                    chrom_length += len(line.strip())

        with open(output.red_ref, 'w') as dump:
            _ = dump.write(out_buffer.getvalue())

    # end of rule


rule build_male_t2t_assembly:
    input:
        t2t = 'references/downloads/T2Tv1_T2TC_chm13.fa.gz',
        hg38 = 'references/assemblies/hg38_GCA_p13.fasta',
    output:
        'references/assemblies/T2Tv1_38p13Y_chm13.fasta'
    resources:
        mem_total_mb = 8192,
        mem_per_cpu_mb = 8192
    run:
        import io
        import gzip

        chrY_buffer = ''
        t2t_buffer = io.StringIO()
        t2t_line_length = 0

        buffer = False
        with open(input.hg38, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    if line.strip() == '>chrY':
                        buffer = True
                        continue
                    elif buffer:
                        # chrY sequence already recorded
                        break
                    else:
                        buffer = False
                        continue
                if buffer:
                    chrY_buffer += line.strip()
        
        chars_written = 0
        with gzip.open(input.t2t, 'rt') as fasta:
            for line in fasta:
                if line.startswith('>chrM'):
                    assert t2t_line_length > 0, 'No T2T line length recorded'
                    t2t_buffer.write('>chrY\n')
                    for i in range(len(chrY_buffer) // t2t_line_length + 1):
                        start = i * t2t_line_length
                        end = start + t2t_line_length
                        chars_written += t2t_buffer.write(chrY_buffer[start:end])
                        _ = t2t_buffer.write('\n')
                    if not chars_written == len(chrY_buffer):
                        raise ValueError('chrY sequence lost: {} vs {}'.format(chars_written, len(chrY_buffer)))
                    t2t_buffer.write(line)
                    continue
                t2t_buffer.write(line)
                if t2t_line_length == 0:
                    if line.strip() and not line.startswith('>'):
                        t2t_line_length = len(line.strip())

        with open(output[0], 'w') as dump:
            _ = dump.write(t2t_buffer.getvalue())
    # end of rule
                    