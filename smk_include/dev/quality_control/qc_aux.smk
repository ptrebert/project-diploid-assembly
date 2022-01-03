
import hashlib

# Define some constant mappings

RS_HIFIEC = 'HASv0161'
RS_HIFIAF = 'ADFv1'
RS_ONTUL = 'GPYv5011SUP'
RS_ONTUL_FILE_EXT = 'guppy-5.0.11-sup-prom_fastq_pass.fastq.gz'
RS_ONTEC = None  # see below

HIFIASM_TIG_TO_KEY = {
    'r_utg': 'UTGRAW',
    'p_utg': 'UTGPRI',
    'p_ctg': 'CTGPRI',
    'a_ctg': 'CTGALT'
}

HIFIASM_KEY_TO_TIG = dict((v, k) for k, v in HIFIASM_TIG_TO_KEY.items())

# derived readsets

def compute_ontec_readset_hash():

    mbg_kmer = config['mbg_init_kmer']
    mbg_window = config['mbg_window_size']
    mbg_resolve = config['mbg_resolve_kmer']
    mbg_combination = f'MBG-k{mbg_kmer}-w{mbg_window}-r{mbg_resolve}'

    # ONTec reads are only created from error-corrected HiFi reads, i.e.,
    # can fix readset information here
    ontec_full_spec = RS_ONTUL + '_' + RS_HIFIEC + '_' + mbg_combination
    ontec_spec_hash = hashlib.md5(ontec_full_spec.encode('ascii')).hexdigest()
    ontec_spec_key = ontec_spec_hash.upper()[:8]
    with open(ontec_spec_key + '.hash.key', 'w') as dump:
        _ = dump.write(ontec_full_spec + '\n')
        _ = dump.write(ontec_spec_hash + '\n')
    return ontec_spec_key


RS_ONTEC = compute_ontec_readset_hash()


def find_script_path(script_name, subfolder=''):
    import os

    current_root = workflow.basedir
    last_root = ''

    script_path = None

    for _ in range(workflow.basedir.count('/')):
        if last_root.endswith('project-diploid-assembly'):
            raise RuntimeError('Leaving project directory tree (next: {}). '
                               'Cannot find script {} (subfolder: {}).'.format(current_root, script_name, subfolder))
        check_path = os.path.join(current_root, 'scripts', subfolder).rstrip('/')  # if subfolder is empty string
        if os.path.isdir(check_path):
            check_script = os.path.join(check_path, script_name)
            if os.path.isfile(check_script):
                script_path = check_script
                break
        last_root = current_root
        current_root = os.path.split(current_root)[0]

    if script_path is None:
        raise RuntimeError('Could not find script {} (subfolder {}). '
                           'Started at path: {}'.format(script_name, subfolder, workflow.basedir))
    return script_path


def validate_readset(readset, input_reads):
    if readset not in input_reads:
        raise ValueError(f'No readset match: {readset} / {input_reads}')
    return None
