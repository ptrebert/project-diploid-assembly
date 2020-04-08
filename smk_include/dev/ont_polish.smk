

rule master:
    input:
        expand('output/helen_polish_clean/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.helpol-p1.fasta',
               hap=['h1', 'h2']),
        expand('output/helen_polish_mixed/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.helmix-p1.fasta',
               hap=['h1', 'h2'])



rule produce_assemblies:
    output:
        protected(expand('input/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.fasta',
                         hap=['h1', 'h2'],
                         num=list(range(1, 25)))),
        protected(expand('input/HG00733_hpg_ontpm-ul_1000.{hap}-un.cluster{num}.fastq.gz',
                         hap=['h1', 'h2'],
                         num=list(range(1, 25))))

rule produce_models:
    output:
        protected('output/helen_models/HELEN_r941_guppy344_human.pkl'),
        protected('output/helen_models/MP_r941_guppy344_human.json')


rule run_alignments:
    input:
        assembly = ancient('input/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.fasta'),
        reads = ancient('input/HG00733_hpg_ontpm-ul_1000.{hap}-un.cluster{num}.fastq.gz')
    output:
        protected('output/alignments/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.bam')
    threads: 48
    shell:
        'minimap2 -ax map-ont -t {threads} {input.assembly} {input.reads} | '
            'samtools sort -@ {threads} | '
            'samtools view -hb -F 0x104 > {output}'


rule samtools_index:
    input:
        '{file_path}.bam'
    output:
        '{file_path}.bam.bai'
    threads: 6
    shell:
        'samtools index -@ {threads} {input}'


rule generate_images:
    input:
        assembly = ancient('input/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.fasta'),
        alignments = ancient('output/alignments/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.bam'),
        bam_index = ancient('output/alignments/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.bam.bai'),
        mp_model = ancient('output/helen_models/MP_r941_guppy344_human.json')
    output:
        directory('output/mp_images/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}/')
    log:
       'log/output/mp_images/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.log'
    threads: 24
    shell:
        'set +u ; '
        'source /home/pebert/work/code/github/external/helen/venv/bin/activate && '
        'marginpolish {input.alignments} {input.assembly} {input.mp_model} '
        '--threads {threads} --produceFeatures --outputBase {output} &> {log} && '
        'deactivate'


rule polish_contigs:
    input:
         img_dir = rules.generate_images.output,
         helen_model = ancient('output/helen_models/HELEN_r941_guppy344_human.pkl')
    output:
          directory('output/helen_polish/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}/')
    log:
       'log/output/helen_polish/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.log'
    threads: 48
    shell:
         'set +u ; '
         'source /home/pebert/work/code/github/external/helen/venv/bin/activate && '
         'helen polish --image_dir {input.img_dir} --model_path {input.helen_model} '
         '--threads {threads} --output_dir {output} --output_prefix HELEN &> {log} && '
         'deactivate'


rule clean_helen_output:
    input:
        'output/helen_polish/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}/'
    output:
        'output/helen_polish_clean/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.fasta',
        'output/helen_polish_clean/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.contigs'
    params:
        contigs = lambda wildcards, input: os.path.join(input[0], 'HELEN.fa')
    run:
        import io

        buffer = io.StringIO()
        polished_contigs = set()

        # HELEN produces single line FASTA output
        with open(params.contigs, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    sequence = line.strip()
                    header = '_'.join([line.strip(), 'cluster' + wildcards.num, wildcards.hap, 'polished'])
                elif not line.strip():
                    continue
                else:
                    if len(line.strip()) > 0:
                        buffer.write(header + '\n' + line)
                        polished_contigs.add(sequence)

        with open(output[0], 'w') as dump:
            _ = dump.write(buffer.getvalue())

        with open(output[1], 'w') as dump:
            _ = dump.write('\n'.join(sorted(polished_contigs)))

rule merge_clean_output:
    input:
        expand('output/helen_polish_clean/HG00733_hpg_ontpm-ul_1000-flye.{{hap}}-un.cluster{num}.fasta',
               num=list(range(1, 25)))
    output:
        'output/helen_polish_clean/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.helpol-p1.fasta',
    shell:
        'cat {input} > {output}'

rule create_mixed_output:
    input:
        raw = 'input/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.fasta',
        clean = 'output/helen_polish_clean/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.fasta',
        polished_contigs = 'output/helen_polish_clean/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.contigs',
    output:
        'output/helen_polish_merged/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.cluster{num}.fasta'
    run:
        import io

        with open(input.polished_contigs, 'r') as known:
            known_contigs = set(known.read().split())

        buffer = io.StringIO()

        with open(input.clean, 'r') as fasta:
            _ = buffer.write(fasta.read())

        with open(input.raw, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    if line.strip() in known_contigs:
                        drop = False
                        line = '_'.join([line.strip(), 'cluster' + wildcards.num, wildcards.hap, 'unpolished']) + '\n'
                    else:
                        drop = True
                    if drop:
                        continue
                    else:
                        buffer.write(line)

        with open(output[0], 'w') as fasta:
            _ = fasta.write(buffer.getvalue())

rule merge_mixed_output:
    input:
         expand('output/helen_polish_merged/HG00733_hpg_ontpm-ul_1000-flye.{{hap}}-un.cluster{num}.fasta',
                num=list(range(1, 25)))
    output:
          'output/helen_polish_mixed/HG00733_hpg_ontpm-ul_1000-flye.{hap}-un.helmix-p1.fasta',
    shell:
         'cat {input} > {output}'