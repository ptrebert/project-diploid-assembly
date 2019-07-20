
include: 'handle_reference_download.smk'

localrules: master_preprocess_references

rule master_preprocess_references:
    input:
        rules.master_handle_reference_download.input,

        expand('references/assemblies/{ref_genome}.fasta',
                ref_genome=config['use_ref_genome'])


rule normalize_reference_assembly_names:
    input:
        'references/downloads/{ref_genome}.fa.gz'
    output:
        seq = 'references/assemblies/{ref_genome}.fasta',
        table = 'references/assemblies/{ref_genome}.sizes'
    run:
        import gzip as gzip
        import io as io

        out_buffer = io.StringIO()
        chromosome_sizes = []
        chrom_length = 0
        current_chrom = None
        with gzip.open(input[0], 'rt') as genome:
            for line in genome:
                if line.startswith('>'):
                    if current_chrom is not None:
                        chromosome_sizes.append((current_chrom, chrom_length))
                        chrom_length = 0
                    chrom_name = line.split()[0].strip('>')
                    try:
                        int(chrom_name)
                        chrom_name = 'chr' + chrom_name
                    except ValueError:
                        if chrom_name in ['X', 'Y', 'M', 'MT']:
                            chrom_name = 'chr' + chrom_name
                    current_chrom = chrom_name
                    out_buffer.write('>{}\n'.format(chrom_name))
                else:
                    out_buffer.write(line)
                    chrom_length += len(line.strip())

        with open(output.seq, 'w') as dump:
            _ = dump.write(out_buffer.getvalue())

        with open(output.table, 'w') as dump:
            for chrom, size in chromosome_sizes:
                _ = dump.write('{}\t{}\n'.format(chrom, size))
    # end of rule