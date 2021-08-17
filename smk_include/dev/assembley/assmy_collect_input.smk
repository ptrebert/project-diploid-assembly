
import pathlib

def read_sample_table(tsv_path):
    import pandas as pd
    import collections as col

    df = pd.read_csv(tsv_path, sep='\t', header=0)
    sample_infos = col.defaultdict(dict)

    for row in df.itertuples(index=False):
        sex = row.sex[0].upper()
        assert sex in ['M', 'F']
        long_sample = f'{row.super_population}-{row.population}-{row.family_id}-{sex}_{row.sample}'
        sample_infos[row.sample]['long_id'] = long_sample
        sample_infos[row.sample]['sex'] = sex
    return sample_infos


def add_hifiec_readsets(sample_infos, hifiec_path):

    fasta_files = pathlib.Path(hifiec_path).glob('*.ec-reads.fasta.gz')
    hifiec_samples = []

    for fasta_file in fasta_files:
        file_size_bytes = fasta_file.stat().st_size
        if file_size_bytes < 15 * 1024**3:
            # assume file is still being dumped to disk
            continue
        file_sample = fasta_file.name.split('_')[0]
        assert file_sample in sample_infos
        sample_infos[file_sample]['HIFIEC'] = fasta_file
        hifiec_samples.append(file_sample)
    return sample_infos, sorted(hifiec_samples)


def add_assembly_graphs(sample_infos, assembly_path):
    import collections as col

    graph_files = pathlib.Path(assembly_path).glob('**/*.gfa')
    assembled_samples = set()
    graph_count = col.Counter()

    for graph_file in graph_files:
        if 'noseq' in graph_file.name:
            continue
        if 'p_utg' in graph_file.name:
            continue
        file_sample = graph_file.name.split('_')[0]
        assert file_sample in sample_infos
        if 'a_ctg' in graph_file.name:
            sample_infos[file_sample]['TIGALT'] = graph_file
        elif 'p_ctg' in graph_file.name:
            sample_infos[file_sample]['TIGPRI'] = graph_file
        elif 'r_utg' in graph_file.name:
            sample_infos[file_sample]['TIGRAW'] = graph_file
        else:
            raise ValueError(f'Unexpected file: {graph_file.name}')
        assembled_samples.add(file_sample)
        graph_count[file_sample] += 1
    for n, c in graph_count.most_common():
        if c != 3:
            raise ValueError(f'Missing assembly graph type for sample {n}')
    return sample_infos, sorted(assembled_samples)


def add_ontul_readsets(sample_infos, ontul_path):

    suffix = 'ONTUL_guppy-5.0.11-sup-prom.fasta.gz'
    fasta_files = pathlib.Path(ontul_path).glob(f'*{suffix}')
    ontul_samples = []

    for fasta_file in fasta_files:
        file_size_bytes = fasta_file.stat().st_size
        if file_size_bytes < 30 * 1024**3:
            # assume file is still incomplete
            continue
        file_sample = fasta_file.name.split('_')[0]
        assert file_sample in sample_infos
        sample_infos[file_sample]['ONTUL'] = fasta_file
        ontul_samples.append(file_sample)

    return sample_infos, sorted(ontul_samples)


def add_ontec_readsets(sample_infos, ontec_path):

    suffix = 'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k1001-w500.clip-ec.fa.gz'
    fasta_files = pathlib.Path(ontec_path).glob(f'*{suffix}')
    ontec_samples = []

    for fasta_file in fasta_files:
        file_size_bytes = fasta_file.stat().st_size
        if file_size_bytes < 30 * 1024**3:
            # assume file is still being dumped to disk
            continue
        file_sample = fasta_file.name.split('_')[0]
        assert file_sample in sample_infos
        sample_infos[file_sample]['ONTEC'] = fasta_file
        ontec_samples.append(file_sample)

    return sample_infos, sorted(ontec_samples)


PATH_SAMPLE_TABLE = config['path_sample_table']
PATH_HIFIEC_READS = config['path_hifiec_reads']
PATH_ASSEMBLY_GRAPHS = config['path_assembly_graphs']
PATH_ONTUL_READS = config['path_ontul_reads']
PATH_ONTEC_READS = config['path_ontec_reads']


def init_samples_and_data():

    location_smk_file = pathlib.Path(workflow.basedir)
    location_sample_table = (location_smk_file / pathlib.Path(PATH_SAMPLE_TABLE)).resolve()
    
    sample_infos = read_sample_table(location_sample_table)
    sample_infos, hifiec_samples = add_hifiec_readsets(sample_infos, PATH_HIFIEC_READS)
    sample_infos, assembled_samples = add_assembly_graphs(sample_infos, PATH_ASSEMBLY_GRAPHS)
    sample_infos, ontul_samples = add_ontul_readsets(sample_infos, PATH_ONTUL_READS)
    sample_infos, ontec_samples = add_ontec_readsets(sample_infos, PATH_ONTEC_READS)

    return sample_infos, hifiec_samples, ontul_samples, ontec_samples, assembled_samples


SAMPLE_INFOS, HIFIEC_SAMPLES, ONTUL_SAMPLES, ONTEC_SAMPLES, ASSEMBLED_SAMPLES = init_samples_and_data()

CONSTRAINT_SAMPLES = '(' + '|'.join(sorted(SAMPLE_INFOS.keys())) + ')'
