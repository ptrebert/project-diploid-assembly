
import os
import csv

import pandas as pd

SEQ_PLATFORMS = {
    'ccs': 'HiFi',
    'clr': 'CLR',
}

SEQ_PLATFORMS_COLORS = {
    'CLR': '#000000',
    'clr': '#000000',
    'HiFi': '#B22222',
    'ccs': '#B22222'
}

POPULATION_ANNOTATION_FILE = '1kg_hgsvc_colors.csv'
SAMPLE_TABLE_FILE = 'sample_table.tsv'

POPULATION_SORT_PLOTTING = [
    'AFR',
    'AMR',
    'EUR',
    'SAS',
    'EAS'
]

PLOT_PROPERTIES = {
    'fontsize_axis_label': 16,
    'fontsize_axis_ticks': 14,
    'fontsize_legend': 14,
    'bar_width': 0.3,
    'bar_intra_spacing': 0.35,
    'bar_inter_spacing': 0.2,
    'CLR_marker': 'D',
    'HiFi_marker': 'o',
    'legend_marker_size': 10,
    'plot_marker_size': 60,
    'dpi_low_res': 150,
    'dpi_high_res': 1200
}


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


def load_plot_data_files(path, extension, pipeline_version):
    """
    """
    search_path = os.path.join(path, pipeline_version)
    data_files = [f for f in os.listdir(search_path) if f.endswith(extension)]
    if not data_files:
        raise ValueError('No data files underneath path {} matching file extension {}'.format(search_path, extension))
    use_files = []
    for df in data_files:
        full_path = os.path.join(search_path, df)
        if df.startswith('v') and not df.startswith(pipeline_version):
            raise ValueError('Wrong version for data file (must be {}): {} / {}'.format(pipeline_version, df, search_path))
        elif df.startswith(pipeline_version):
            new_name = df.split('_', 1)[1]
            new_path = os.path.join(search_path, new_name)
            os.rename(full_path, new_path)
            use_files.append(new_path)
        else:
            use_files.append(full_path)
    return sorted(use_files)


def extract_sample_platform(filename, multi_readset=False, mapped_readset=False, long_read_pos=2):
    """
    """
    if multi_readset:
        parts = filename.split('.')
        sample = parts[0]
        lr_data = parts[long_read_pos].split('_')
        platform = lr_data[1].split('-')[-1]
    elif mapped_readset:
        parts = filename.split('_map-to_')
        sample, readset1 = parts[0].split('_', 1)
        readset2 = parts[-1].split('.')[0]
        if long_read_pos == 1:
            lr_data = readset1.split('_')
        elif long_read_pos == 2:
            lr_data = readset2.split('_')
        else:
            raise ValueError('For aligned readsets, postion of long read info has to be 1 or 2, but not {}'.format(long_read_pos))
        platform = lr_data[1].split('-')[-1]
    else:
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


def load_sample_dataframe(tsv_path=None):
    """
    """
    if tsv_path is None or not os.path.isfile(tsv_path):
        tsv_path = find_annotation_file(SAMPLE_TABLE_FILE)
    df = pd.read_csv(tsv_path, sep='\t')
    return df


def check_cache_consistency(plot_data):
    """
    """
    samples =  load_sample_dataframe()
    skip_samples = samples.loc[samples['2020_SKIP'] == 1, 'individual']
    skip_samples = sorted([s for s in skip_samples if s in plot_data.index.get_level_values('sample')])

    missing_samples = []
    use_samples = samples['2020_SKIP'] == 0
    for platform in SEQ_PLATFORMS.values():
        plot_subset = plot_data.xs(platform, level='platform')
        subset_samples = set(plot_subset.index.get_level_values('sample'))
        
        platform_samples = samples[platform] == 1
        expected_samples = samples.loc[use_samples & platform_samples, 'individual']
        [missing_samples.append((s, platform)) for s in expected_samples if s not in subset_samples]
    
    return skip_samples, missing_samples


def get_grey_bg(hex=False):
    if hex:
        grey_color = '#b4b4b4'.upper()
    else:
        rgb = 180, 180, 180
        rgb = map(lambda x: round(x/255, 2), rgb)
        grey_color = tuple(rgb)
    return grey_color


def hex_to_rgb(hex_string, norm=True):

    h = hex_string.strip('#')
    if norm:
        rgb = tuple(x/255 for x in [int(h[i:i+2], 16) for i in [0, 2, 4]])
    else:
        rgb = tuple(int(h[i:i+2], 16) for i in [0, 2, 4])
    return rgb


def relative_rgb(rgb):
    return tuple(round(x/255, 2) for x in rgb)


def get_gray_bg(hex=False):
    return get_grey_bg(hex)


def get_platform_color(platform, hex=False):
    assert platform in SEQ_PLATFORMS or platform in SEQ_PLATFORMS.values(), 'Unknown platform: {}'.format(platform)
    if hex:
        return SEQ_PLATFORMS_COLORS[platform]
    else:
        return hex_to_rgb(SEQ_PLATFORMS_COLORS[platform])


def get_population_sorting():
    return POPULATION_SORT_PLOTTING


def get_plot_property(key):
    try:
        prop = PLOT_PROPERTIES[key]
    except KeyError:
        known_keys = sorted(PLOT_PROPERTIES.keys())
        raise KeyError('Unknown key: {} (known keys: {})'.format(key, known_keys))
    return prop


def get_sequencing_platforms():
    """
    """
    return sorted(SEQ_PLATFORMS.values())


def add_incomplete_stamp(mpl_axis, xloc, yloc):
    """
    """
    mpl_axis.text(
            x=xloc,
            y=yloc,
            s='incomplete',
            color='red',
            rotation=-45,
            alpha=0.5,
            fontdict={
                'weight': 'bold',
                'size': 16
            },
            transform=mpl_axis.transAxes,
            zorder=3
        )
    return mpl_axis
