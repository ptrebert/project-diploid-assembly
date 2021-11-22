localrules: dump_mbg_param_info

import hashlib


HIFIASM_TIGS = {
    'UTGRAW': 'r_utg',
    'UTGPRI': 'p_utg',
    'CTGPRI': 'p_ctg',
    'CTGALT': 'a_ctg'
}

def get_hifiasm_tigs(tig_spec):

    try:
        result = HIFIASM_TIGS[tig_spec]
    except KeyError:
        try:
            result = HIFIASM_TIGS[tig_spec.split('-')[1]]
        except KeyError:
            raise ValueError(f'Unknown MBG param spec: {tig_spec} / {HIFIASM_TIGS}')
    return result


def compute_lut_mbg_params():

    init_kmer = config['mbg_init_kmer']
    window_sizes = config['mbg_window_size']
    resolve_kmer = config['mbg_resolve_kmer']

    if not len(init_kmer) == len(window_sizes) == len(resolve_kmer):
        raise ValueError('The same number of parameters has to be specified for MBG k/w/r')

    param_lut = dict()
    for k, w, r in zip(init_kmer, window_sizes, resolve_kmer):
        param_comb = f'k{k}-w{w}-r{r}'.encode('ascii')
        param_hash = hashlib.sha1(param_comb).hexdigest()
        key = param_hash.upper()[:6]
        if key in param_lut:
            raise ValueError(f"MBG param key collision: {param_comb.decode('ascii')} / {param_lut[key]}")
        param_lut[key] = k, w, r
    return param_lut


MBG_PARAMS = compute_lut_mbg_params()


def get_mbg_param(param_spec, which=None):

    param_pos = {
        'k': 0,
        'w': 1,
        'r': 2,
        'kmer': 0,
        'window': 1,
        'resolve': 2
    }

    try:
        result = MBG_PARAMS[param_spec]
    except KeyError:
        try:
            result = MBG_PARAMS[param_spec.split('-')[1]]
        except KeyError:
            raise ValueError(f'Unknown MBG param spec: {param_spec} / {MBG_PARAMS}')
    if which is not None:
        assert which in param_pos, f'Cannot process "which": {which}'
        result = result[param_pos[which]]
    return result


rule dump_mbg_param_info:
    output:
        expand(
            'assembler_params/MBG_{param_info}.info',
            param_info=['{}_k{}-w{}-r{}'.format(phash, *pvalues) for phash, pvalues in MBG_PARAMS.items()]
        )
    priority: 100
    run:
        file_template = 'assembler_params/MBG_{}.info'
        for phash, pvalues in MBG_PARAMS.items():
            param_info = f'{phash}_{pvalues}'
            with open(file_template.format(param_info), 'w') as dump:
                pass


def compute_lut_lja_params():

    small_kmer = config['lja_small_kmer']
    large_kmer = config['lja_large_kmer']

    if not len(small_kmer) == len(large_kmer):
        raise ValueError('The same number of parameters has to be specified for LJA k/K')

    param_lut = dict()
    for k, K in zip(small_kmer, large_kmer):
        param_comb = f'k{k}-K{K}'.encode('ascii')
        param_hash = hashlib.sha1(param_comb).hexdigest()
        key = param_hash.upper()[:6]
        if key in param_lut:
            raise ValueError(f"LJA param key collision: {param_comb.decode('ascii')} / {param_lut[key]}")
        param_lut[key] = k, K
    return param_lut


LJA_PARAMS = compute_lut_lja_params()


def get_lja_param(param_spec, which=None):

    param_pos = {
        'k': 0,
        'K': 1,
        'smallk': 0,
        'largek': 1
    }

    try:
        result = LJA_PARAMS[param_spec]
    except KeyError:
        try:
            result = LJA_PARAMS[param_spec.split('-')[1]]
        except KeyError:
            raise ValueError(f'Unknown LJA param spec: {param_spec} / {LJA_PARAMS}')
    if which is not None:
        assert which in param_pos, f'Cannot process "which": {which}'
        result = result[param_pos[which]]
    return result