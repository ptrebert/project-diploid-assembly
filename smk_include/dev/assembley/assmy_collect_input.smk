localrules: link_input_hifi_ec_reads,
            link_input_hifi_assemblies


ROOT_FOLDER_SAMPLE_CONFIGS = '/home/local/work/code/github/project-diploid-assembly/smk_config/samples/hifi_v13'
ROOT_FOLDER_ASSM_GRAPH = ''
ROOT_FOLDER_HIFI_EC_READS = ''


def collect_sexed_samples(top_dir):
    import os
    import yaml
    import pathlib

    males = set()
    females = set()
    sample_long_id = dict()

    for root, dirs, files in os.walk(top_dir):
        sample_configs = [f for f in files if f.endswith('.yml')]
        for s in sample_configs:
            sample = s.split('.')[0].upper()
            with open(pathlib.Path(root, s), 'rb') as yaml_dump:
                cfg = yaml.load(yaml_dump)
                sample_info = cfg['sample_description_' + sample]
                assert sample_info['individual'] == sample, f'Sample config error: {sample} / {s} / {sample_info}'
                sex = sample_info['sex']
                pop = sample_info['ASK']
                super_pop = sample_info['super_population']
                long_sample = f'{super_pop}_{pop}_{sample}'
                if sex == 'male':
                    males.add(sample)
                    sample_long_id[sample] = long_sample
                elif sex == 'female':
                    females.add(sample)
                    sample_long_id[sample] = long_sample
                else:
                    raise ValueError(f'Unknown sex: {sample_info}')

    return sorted(males), sorted(females), sample_long_id


def find_hifi_ec_reads(input_dir, samples):
    import os
    import pathlib

    read_files = os.listdir(input_dir)
    read_files = [f for f in read_files if f.split('_')[0] in samples]
    read_files = [str(pathlib.Path(input_dir, f)) for f in read_files]

    return sorted(read_files)


rule link_input_hifi_ec_reads:
    input:
        reads = find_hifi_ec_reads
    output:
        link = 'input/ec_reads/{sample}.hifi.fasta.gz'
    run:
        import os
        short_sample = wildcards.sample.split('_')[-1]

        for read_file in input.reads:
            if not read_file.startswith(short_sample):
                continue
            os.symlink(read_file, output.link)
            break




males, females, sample_long_id = collect_sexed_samples(ROOT_FOLDER_SAMPLE_CONFIGS)
