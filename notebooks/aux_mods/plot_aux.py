
import os
import csv

SEQ_PLATFORMS = {
    'ccs': 'HiFi',
    'clr': 'CLR',
}

POPULATION_ANNOTATION_FILE = '1kg_hgsvc_colors.csv'
SAMPLE_TABLE_FILE = 'sample_table.tsv'

POPULATION_SORT_PLOTTING = [
    'AFR',
    'SAS',
    'EAS',
    'EUR',
    'AMR'
]


def find_annotation_file(file_name, search_path=None):
    """
    """
    if search_path is None:
        search_path = __name__
    search_path = os.path.abspath(search_path)

    source_path = search_path
    sample_table_path = ''
    while 1:
        check_path = os.path.join(search_path, 'annotation', file_name)
        if os.path.isfile(check_path):
            sample_table_path = check_path
            break
        search_path = os.path.split(search_path)[0]
        if search_path == '/' or not search_path:
            raise ValueError('Could not find sample table starting at: {}'.format(source_path))
    return sample_table_path


def load_population_annotation(tsv_path=None, relative_rgb=True):
    """
    """
    if tsv_path is None or not os.path.isfile(tsv_path):
        tsv_path = find_annotation_file(POPULATION_ANNOTATION_FILE)
    pop_to_hex = {}
    pop_to_rgb = {}
    pop_to_super = {}
    with open(tsv_path, 'r', newline='') as table:
        reader = csv.DictReader(table, delimiter='\t')
        for record in reader:
            pop_to_hex[record['population']] = '#' + record['hex'].upper()
            rgb = tuple(map(int, record['rgb'].split(',')))
            if relative_rgb:
                rgb = tuple(round(x/255, 2) for x in rgb)
            pop_to_rgb[record['population']] = rgb
            pop_to_super[record['population']] = record['super_pop']
    return pop_to_hex, pop_to_rgb, pop_to_super


def extract_sample_platform(filename):
    """
    """
    parts = filename.split('_')
    sample = parts[0]
    platform = parts[2].split('-')[-1]
    return sample, SEQ_PLATFORMS[platform]


def load_sample_table(tsv_path=None):
    """
    """
    if tsv_path is None or not os.path.isfile(tsv_path):
        tsv_path = find_annotation_file(SAMPLE_TABLE_FILE)
    sample_lut = {}
    with open(tsv_path, 'r', newline='') as table:
        reader = csv.DictReader(table, delimiter='\t')
        for record in reader:
            sample_lut[record['individual']] = record
    return sample_lut


def get_grey_bg(hex=False):
    if hex:
        grey_color = '#b4b4b4'.upper()
    else:
        rgb = 180, 180, 180
        rgb = map(lambda x: round(x/255, 2), rgb)
        grey_color = tuple(rgb)
    return grey_color


def get_gray_bg(hex=False):
    return get_grey_bg(hex)


def get_population_sorting():
    return POPULATION_SORT_PLOTTING
