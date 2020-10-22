
localrules: master_bng_hybrids

HYBRID_CONFIG = {
    'chroms': ['chr' + str(i) for i in range(1, 23)] + ['chrXY', 'chrUn'],
    'reference': 'GRCh38_HGSVC2_noalt'
}


def bng_hybrids_determine_targets(wildcards):

    hybrid_targets = {
        'contig_stats': 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.contig-stats.tsv',
        'scaffold_align': 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}.bed',
        'unsupport_align': 'output/evaluation/bng_hybrids/{assembly}/{assembly}.unsupported_map-to_{reference}.bed',
        'merged_align': 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid_map-to_{reference}.merged.bed'
    }

    search_paths = {
        (os.path.join(os.getcwd(), 'output/evaluation/scaffolded_assemblies'), 'contig_stats'),
        (os.path.join(os.getcwd(), 'output/evaluation/scaffolded_assemblies'), 'scaffold_align'),
        (os.path.join(os.getcwd(), 'output/evaluation/scaffolded_assemblies'), 'unsupport_align'),
        (os.path.join(os.getcwd(), 'output/evaluation/scaffolded_assemblies'), 'merged_align')
    }

    fixed_wildcards = {
        'reference': HYBRID_CONFIG['reference']
    }

    compute_results = set()

    for path, trg_type in search_paths:
        if not os.path.isdir(path):
            continue
        hybrid_trg = hybrid_targets[trg_type]
        for hybrid_file in os.listdir(path):
            if not hybrid_file.endswith('.bng-hybrid.agp'):
                continue
            tmp = dict(fixed_wildcards)
            if trg_type in ['unsupport_align', 'merged_align']:
                tmp['assembly'] = hybrid_file.rsplit('.', 2)[0]
                fmt_target = hybrid_trg.format(**tmp)
                compute_results.add(fmt_target)
            else:
                for c in HYBRID_CONFIG['chroms']:
                    tmp['assembly'] = hybrid_file.rsplit('.', 2)[0]
                    tmp['chrom'] = c
                    fmt_target = hybrid_trg.format(**tmp)
                    compute_results.add(fmt_target)

    return sorted(compute_results)


rule master_bng_hybrids:
    input:
        bng_hybrids_determine_targets


rule summarize_hybrid_statistics:
    input:
        agp = 'output/evaluation/scaffolded_assemblies/{assembly}.bng-hybrid.agp',
        fasta = 'output/evaluation/scaffolded_assemblies/{assembly}.bng-scaffolds.fasta',
        discard = 'output/evaluation/scaffolded_assemblies/{assembly}.bng-unsupported.fasta',
        ctg_ref_aln = 'output/alignments/contigs_to_reference/evaluation/phased_assemblies/{assembly}_map-to_GRCh38_HGSVC2_noalt.bed',
        dummy_fasta = 'references/assemblies/GRCh38_HGSVC2_noalt.dummy.fasta',
        assm_index = 'output/evaluation/phased_assemblies/{assembly}.fasta.fai'
    output:
        contig_stats = 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.contig-stats.tsv',
        layout = 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.scaffold-layout.tsv',
        scaffolds = 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.scaffolds.wg.fasta',
        contig_fasta = expand(
            'output/evaluation/bng_hybrids/{{assembly}}/{{assembly}}.hybrid.contigs.{chrom}.fasta',
            chrom=HYBRID_CONFIG['chroms']
        )
    log:
        'log/output/evaluation/bng_hybrids/summarize/{assembly}.hybrid.stats.log',
    conda:
        '../../environment/conda/conda_evaltools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 12288 * attempt,
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
    params:
        exec = lambda wildcards: find_script_path('process_bng_hybrid.py'),
        out_prefix = lambda wildcards, output: output.contig_stats.rsplit('.', 3)[0] + '.hybrid'
    shell:
        '{params.exec} --agp-file {input.agp} --fasta-file {input.fasta} --dummy-fasta {input.dummy_fasta} '
            '--bed-file {input.ctg_ref_aln} --assembly-index {input.assm_index} --output {params.out_prefix} &> {log}'


def write_fasta(entry, sequence, output_file):

    import io
    out_buffer = io.StringIO()

    line_width = 120
    lines = len(sequence) // line_width + 1
    chars_written = 0

    _ = out_buffer.write('>{}\n'.format(entry))

    for pos in range(lines):
        chars_written += out_buffer.write(sequence[pos * line_width:pos * line_width + line_width])
        _ = out_buffer.write('\n')

    _ = out_buffer.write('\n')

    if not chars_written == len(sequence):
        raise ValueError('Lost parts of sequence: {} / {} / {} / {} / {}'.format(entry, len(sequence), chars_written, line_width, lines))

    _ = output_file.write(out_buffer.getvalue())
    return


rule split_reference_assembly:
    input:
        'references/assemblies/{reference}.fasta'
    output:
        expand(
            'references/assemblies/{{reference}}.{seq_parts}.fasta',
            seq_parts=['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrXY', 'wg-male', 'wg-female']
        ),
        'references/assemblies/{reference}.dummy.fasta',

    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 12288 * attempt,
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt
    run:
        autosomes = ['chr' + str(i) for i in range(1, 23)]
        female = ['chrX']
        male = ['chrY']
        primary_chroms = autosomes + female + male

        current_chrom = ''
        current_seq = ''
        seq_buffer = dict()
        buffer = False

        with open(input[0], 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    chrom = line.strip().strip('>')
                    if chrom not in primary_chroms:
                        if current_chrom and buffer:
                            seq_buffer[current_chrom] = current_seq
                            current_chrom = ''
                            current_seq = ''
                        buffer = False
                        continue
                    if current_chrom and buffer:
                        seq_buffer[current_chrom] = current_seq
                    current_chrom = chrom
                    current_seq = ''
                    buffer = True
                    continue
                current_seq += line.strip()
        
        if buffer:
            seq_buffer[current_chrom] = current_chrom

        for c in autosomes + female:
            outfile = input[0].replace('.fasta', '.{}.fasta')
            tmp = outfile.format(c)
            with open(tmp, 'w') as dump:
                write_fasta(c, seq_buffer[c], dump)

        male_chroms = input[0].replace('.fasta', '.chrXY.fasta')
        with open(male_chroms, 'w') as dump:
            write_fasta('chrX', seq_buffer['chrX'], dump)
            write_fasta('chrY', seq_buffer['chrY'], dump)

        male_genome = input[0].replace('.fasta', '.wg-male.fasta')
        with open(male_genome, 'w') as dump:
            for c in autosomes + female + male:
                write_fasta(c, seq_buffer[c], dump)

        female_genome = input[0].replace('.fasta', '.wg-female.fasta')
        with open(female_genome, 'w') as dump:
            for c in autosomes + female:
                write_fasta(c, seq_buffer[c], dump)
        
        # create one file containing a dummy sequence that will
        # be used to ensure non-empty outputs to make Snakemake happy
        dummy = input[0].replace('.fasta', '.dummy.fasta')
        with open(dummy, 'w') as dump:
            write_fasta('dummy', seq_buffer['chr1'][1000000:1500000], dump)
    ### END OF RUN BLOCK


def select_human_reference(wildcards):

    ref_genome = ''
    use_chrom = None

    # determine sample sex
    sample = wildcards.assembly.split('_')[0]
    sample_cfg = config['sample_description_{}'.format(sample)]
    sample_sex = sample_cfg['sex']

    if hasattr(wildcards, 'chrom') and wildcards.chrom.lower() in ['chrun', 'chrxy']:
        if sample_sex == 'male':
            if wildcards.chrom.lower() == 'chrxy':
                use_chrom = 'chrXY'
            else:
                use_chrom = 'wg-male'
        elif sample_sex == 'female':
            if wildcards.chrom.lower() == 'chrxy':
                use_chrom = 'chrX'
            else:
                use_chrom = 'wg-female'
        else:
            raise ValueError('Unrecognized sex: {} / {}'.format(sample_sex, sample))
    elif not hasattr(wildcards, 'chrom'):
        if sample_sex == 'male':
            use_chrom = 'wg-male'
        elif sample_sex == 'female':
            use_chrom = 'wg-female'
        else:
            raise ValueError('Unrecognized sex: {} / {}'.format(sample_sex, sample))
    else:
        use_chrom = wildcards.chrom    
    ref_genome = 'references/assemblies/{{reference}}.{}.fasta'.format(use_chrom)

    return ref_genome


rule minimap_scaffold_to_reference_alignment:
    input:
        contigs = 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.contigs.{chrom}.fasta',
        reference = select_human_reference
    output:
        'output/alignments/scaffolds_to_reference/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}.psort.sam.bam'
    log:
        minimap = 'log/output/alignments/scaffolds_to_reference/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}.minimap.log',
        st_sort = 'log/output/alignments/scaffolds_to_reference/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}.st-sort.log',
        st_view = 'log/output/alignments/scaffolds_to_reference/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}.st-view.log',
    benchmark:
        '.'.join(['run/output/alignments/scaffolds_to_reference/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}', 't{}'.format(config['num_cpu_high']), 'rsrc'])
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((8192 + 8192 * attempt) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 8192 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_sort_mb = 4096
    params:
        individual = lambda wildcards: wildcards.assembly.split('_')[0],
        readgroup_id = lambda wildcards: wildcards.assembly.replace('.', '_'),
        discard_flag = config['minimap_contigref_aln_discard'],
        tempdir = lambda wildcards: os.path.join(
                                        'temp', 'minimap', wildcards.reference,
                                        wildcards.assembly, wildcards.chrom)
    shell:
        'rm -rfd {params.tempdir} ; mkdir -p {params.tempdir} && '
        'minimap2 -t {threads} '
            '--secondary=no --eqx -Y -ax asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 '
            '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" '
            '{input.reference} {input.contigs} 2> {log.minimap} | '
            'samtools sort -m {resources.mem_sort_mb}M -T {params.tempdir} 2> {log.st_sort} | '
            'samtools view -b -F {params.discard_flag} /dev/stdin > {output} 2> {log.st_view}'


rule dump_scaffold_to_reference_alignment_to_bed:
    input:
        'output/alignments/scaffolds_to_reference/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}.psort.sam.bam'
    output:
        'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}.bed'
    log:
        'log/output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid_map-to_{reference}.{chrom}.dump-bed.log'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        chrom = '(' + '|'.join(HYBRID_CONFIG['chroms']) + ')'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_total_mb = 2048,
        mem_per_cpu_mb = 2048
    shell:
        'bedtools bamtobed -i {input} > {output} 2> {log}'


def split_contig_to_reads(ctg_name, ctg_seq, buffer):

    if len(ctg_seq) < 1001:
        # in case there are a lot of short contigs,
        # avoid creating too many (overly) short reads
        midpoint = len(ctg_seq) // 2
        seq_a, seq_b = ctg_seq[:midpoint], ctg_seq[midpoint:]
        buffer.write('>{}_read0\n'.format(ctg_name))
        buffer.write('{}\n\n'.format(seq_a))
        buffer.write('>{}_read1\n'.format(ctg_name))
        buffer.write('{}\n\n'.format(seq_b))
    else:
        for i in range(len(ctg_seq) // 500):
            buffer.write('>{}_read{}\n'.format(ctg_name, i))
            buffer.write('{}\n\n'.format(ctg_seq[i*500:i*500+500]))
    return buffer


rule split_unsupported_contigs:
    input:
        'output/evaluation/scaffolded_assemblies/{assembly}.bng-unsupported.fasta',
    output:
        'output/evaluation/scaffolded_assemblies/ctg_reads/{assembly}.unscf-reads.fasta',
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_total_mb = 2048,
        mem_per_cpu_mb = 2048
    run:
        import io

        buffer = io.StringIO()

        current_ctg = None
        current_seq = ''
        with open(input[0], 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    ctg_name = line.strip().strip('>')
                    if current_ctg is not None:
                        if not len(current_seq) < 500:
                            buffer = split_contig_to_reads(current_ctg, current_seq, buffer)
                    current_ctg = ctg_name
                    current_seq = ''
                    continue
                current_seq += line.strip()
        if not len(current_seq) < 500:
            buffer = split_contig_to_reads(current_ctg, current_seq, buffer)

        with open(output[0], 'w') as dump:
            _ = dump.write(buffer.getvalue())
    ### END OF run block


rule minimap_unscaffolded_to_reference_alignment:
    input:
        contigs = 'output/evaluation/scaffolded_assemblies/ctg_reads/{assembly}.unscf-reads.fasta',
        reference = select_human_reference
    output:
        'output/alignments/scaffolds_to_reference/{assembly}/{assembly}.unsupported_map-to_{reference}.psort.sam.bam'
    log:
        minimap = 'log/output/alignments/scaffolds_to_reference/{assembly}/{assembly}.unsupported_map-to_{reference}.minimap.log',
        st_sort = 'log/output/alignments/scaffolds_to_reference/{assembly}/{assembly}.unsupported_map-to_{reference}.st-sort.log',
        st_view = 'log/output/alignments/scaffolds_to_reference/{assembly}/{assembly}.unsupported_map-to_{reference}.st-view.log',
    benchmark:
        '.'.join(['run/output/alignments/scaffolds_to_reference/{assembly}/{assembly}.unsupported_map-to_{reference}', 't{}'.format(config['num_cpu_high']), 'rsrc'])
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((12288 + 12288 * attempt) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 12288 + 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_sort_mb = 4096
    params:
        individual = lambda wildcards: wildcards.assembly.split('_')[0],
        readgroup_id = lambda wildcards: wildcards.assembly.replace('.', '_'),
        discard_flag = config['minimap_contigref_aln_discard'],
        tempdir = lambda wildcards: os.path.join(
                                        'temp', 'minimap', wildcards.reference,
                                        wildcards.assembly)
    shell:
        'rm -rfd {params.tempdir} ; mkdir -p {params.tempdir} && '
        'minimap2 -t {threads} '
            '--secondary=no --eqx -Y -a -x sr '
            '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" '
            '{input.reference} {input.contigs} 2> {log.minimap} | '
            'samtools sort -m {resources.mem_sort_mb}M -T {params.tempdir} 2> {log.st_sort} | '
            'samtools view -b -F {params.discard_flag} /dev/stdin > {output} 2> {log.st_view}'


rule dump_unsupported_to_reference_alignment_to_bed:
    input:
        'output/alignments/scaffolds_to_reference/{assembly}/{assembly}.unsupported_map-to_{reference}.psort.sam.bam'
    output:
        'output/evaluation/bng_hybrids/{assembly}/{assembly}.unsupported_map-to_{reference}.bed'
    log:
        'log/output/evaluation/bng_hybrids/{assembly}/{assembly}.unsupported_map-to_{reference}.dump-bed.log'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_total_mb = 2048,
        mem_per_cpu_mb = 2048
    shell:
        'bedtools bamtobed -i {input} > {output} 2> {log}'


rule merge_scaffold_alignments:
    input:
        expand('output/evaluation/bng_hybrids/{{assembly}}/{{assembly}}.hybrid_map-to_{{reference}}.{chrom}.bed',
                chrom=[c for c in HYBRID_CONFIG['chroms'] if c != 'chrUn'])  # keep chrUn separate to avoid mix-ups
    output:
        'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid_map-to_{reference}.merged.bed'
    shell:
        'cat {input} > {output}'
