
localrules: link_hgsvc_sequel2_ccs_fastq_data, link_hgsvc_sequel2_ccs_bam_data

checkpoint link_hgsvc_sequel2_ccs_fastq_data:
    output:
        HG00731 = expand('input/fastq/partial/parts/HG00731_hgsvc_pbsq2-ccs.part{num}.fastq.gz',
                            num=list(range(1, 6))),
        HG00732 = expand('input/fastq/partial/parts/HG00732_hgsvc_pbsq2-ccs.part{num}.fastq.gz',
                            num=list(range(1, 7))),
        HG00733 = expand('input/fastq/partial/parts/HG00733_hgsvc_pbsq2-ccs.part{num}.fastq.gz',
                            num=list(range(1, 8)))
    run:
        base_input = config['link_base_input']

        hg00731_files = [
            'HG00731_20190925_EEE_m54329U_190528_231241.Q20.fastq.gz',
            'HG00731_20190925_EEE_m54329U_190531_175004.Q20.fastq.gz',
            'HG00731_20190925_EEE_m54329U_190601_235636.Q20.fastq.gz',
            'HG00731_20190925_EEE_m54329U_190606_045526.Q20.fastq.gz',
            'HG00731_20190925_EEE_m54329U_190906_205127.Q20.fastq.gz'
        ]

        hg00732_files = [
            'HG00732_20190925_EEE_m54329U_190604_224858.Q20.fastq.gz',
            'HG00732_20190925_EEE_m54329U_190610_071123.Q20.fastq.gz',
            'HG00732_20190925_EEE_m54329U_190611_132246.Q20.fastq.gz',
            'HG00732_20190925_EEE_m54329U_190612_193411.Q20.fastq.gz',
            'HG00732_20190925_EEE_m54329U_190705_000551.Q20.fastq.gz',
            'HG00732_20190925_EEE_m54329U_190901_130210.Q20.fastq.gz'
        ]

        hg00733_files = [
            'HG00733_20190925_EEE_m54329U_190607_185248.Q20.fastq.gz',
            'HG00733_20190925_EEE_m54329U_190615_010947.Q20.fastq.gz',
            'HG00733_20190925_EEE_m54329U_190617_231905.Q20.fastq.gz',
            'HG00733_20190925_EEE_m54329U_190619_052546.Q20.fastq.gz',
            'HG00733_20190925_EEE_m54329U_190629_180018.Q20.fastq.gz',
            'HG00733_20190925_EEE_m54329U_190701_222759.Q20.fastq.gz',
            'HG00733_20190925_EEE_m54329U_190827_173812.Q20.fastq.gz'
        ]
        import os
        work_dir = os.getcwd()
        link_info = [
            ('HG00731', hg00731_files, sorted(output.HG00731)),
            ('HG00732', hg00732_files, sorted(output.HG00732)),
            ('HG00733', hg00733_files, sorted(output.HG00733))
        ]

        for individual, sources, targets in link_info:
            source_dir = os.path.join(base_input, individual)
            for src, trg in zip(sources, targets):
                source_file = os.path.join(source_dir, src)
                assert os.path.isfile(source_file), 'Path is not a file: {}'.format(source_file)
                target_file = os.path.join(work_dir, trg)
                if not os.path.islink(target_file):
                    os.symlink(source_file, target_file)


checkpoint link_hgsvc_sequel2_ccs_bam_data:
    output:
        HG00731 = expand('input/bam/partial/parts/HG00731_hgsvc_pbsq2-ccs.part{num}.pbn.bam',
                            num=list(range(1, 6))),
        HG00732 = expand('input/bam/partial/parts/HG00732_hgsvc_pbsq2-ccs.part{num}.pbn.bam',
                            num=list(range(1, 7))),
        HG00733 = expand('input/bam/partial/parts/HG00733_hgsvc_pbsq2-ccs.part{num}.pbn.bam',
                            num=list(range(1, 8)))
    run:
        base_input = config['link_base_input']

        hg00731_files = [
            'HG00731_20190925_EEE_m54329U_190528_231241.ccs.bam',
            'HG00731_20190925_EEE_m54329U_190531_175004.ccs.bam',
            'HG00731_20190925_EEE_m54329U_190601_235636.ccs.bam',
            'HG00731_20190925_EEE_m54329U_190606_045526.ccs.bam',
            'HG00731_20190925_EEE_m54329U_190906_205127.ccs.bam'
        ]

        hg00732_files = [
            'HG00732_20190925_EEE_m54329U_190604_224858.ccs.bam',
            'HG00732_20190925_EEE_m54329U_190610_071123.ccs.bam',
            'HG00732_20190925_EEE_m54329U_190611_132246.ccs.bam',
            'HG00732_20190925_EEE_m54329U_190612_193411.ccs.bam',
            'HG00732_20190925_EEE_m54329U_190705_000551.ccs.bam',
            'HG00732_20190925_EEE_m54329U_190901_130210.ccs.bam'
        ]

        hg00733_files = [
            'HG00733_20190925_EEE_m54329U_190607_185248.ccs.bam',
            'HG00733_20190925_EEE_m54329U_190615_010947.ccs.bam',
            'HG00733_20190925_EEE_m54329U_190617_231905.ccs.bam',
            'HG00733_20190925_EEE_m54329U_190619_052546.ccs.bam',
            'HG00733_20190925_EEE_m54329U_190629_180018.ccs.bam',
            'HG00733_20190925_EEE_m54329U_190701_222759.ccs.bam',
            'HG00733_20190925_EEE_m54329U_190827_173812.ccs.bam'
        ]
        import os
        work_dir = os.getcwd()
        link_info = [
            ('HG00731', hg00731_files, sorted(output.HG00731)),
            ('HG00732', hg00732_files, sorted(output.HG00732)),
            ('HG00733', hg00733_files, sorted(output.HG00733))
        ]

        for individual, sources, targets in link_info:
            source_dir = os.path.join(base_input, individual)
            for src, trg in zip(sources, targets):
                source_file = os.path.join(source_dir, src)
                assert os.path.isfile(source_file), 'Path is not a file: {}'.format(source_file)
                target_file = os.path.join(work_dir, trg)
                if not os.path.islink(target_file):
                    os.symlink(source_file, target_file)
