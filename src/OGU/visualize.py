import argparse
import re
from pathlib import Path


import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from pycirclize import Circos
from pycirclize.parser import Genbank

from OGU.global_vars import log
from OGU.gb2fasta import gb2fasta_main

# Evaluation.csv -> xlsx -> first column Type, second column name
part_color = '#6BCB77,#4D96FF,#6BCB77,#FF6B6B'.split(',')
# yellow,yellow/gray,blue,blue,green?
track_colors = list(reversed('#F7931B,#66B3FF,#338CFF,#3F51B5,#FB9883,#CC5500,'
                             '#FF6347,#36454F,#778877'.split(',')))
track_colors2 = ('#F7931B,#66B3FF,#338CFF,#3F51B5,#FB9883,#CC5500,#FF6347,'
                 '#36454F,#778877'.split(','))
colname_text = (('Tree_Res', 'Tree resolution'), ('PD', 'PD'),
                ('PD_stem', 'PD-stem'), ('PD_terminal', 'PD-terminal'),
                ('Observed_Res', 'Observed resolution'), ('Pi', 'Pi'),
                ('Entropy', 'Shannon index'),
                ('Total_GC', 'GC ratio'),
                ('Gap_Ratio', 'Gap ratio'))


class MyRadius:
    def __init__(self):
        self._radius = 100

    def get(self, thin=False):
        width = 4 if thin else 8
        value = self._radius - width, self._radius
        self._radius -= width
        return value


def parse_args(arg_list=None):
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=visualize_main.__doc__)
    arg.add_argument('-input_csv', help='Evaluation result tabel, csv file')
    arg.add_argument('-ref_gb', help='reference genome genbank file')
    arg.add_argument('-taxon', help='Taxa name, recommend use family name')
    arg.add_argument('-type', choices=('cp', 'mt'), dest='og_type',
                     help='Type of organelle (mt:mitochondria or '
                          'cp:plastid/chloroplast')
    arg.add_argument('-n_seqs', type=int, dest='count_threshold', default=10,
                     help='Minimum number of sequences per gene/fragment')
    # use tobacco as default
    arg.add_argument('-lsc', type=int, default=86686,
                     help='LSC size of plastid genome')
    arg.add_argument('-ssc', type=int, default=18571,
                     help='SSC size of plastid genome')
    arg.add_argument('-ir', type=int, default=25341,
                     help='IR size of plastid genome')
    arg.add_argument('-out', help='Output file', default='figure.pdf')
    if arg_list is None:
        return arg.parse_known_args()
    else:
        return arg.parse_known_args(arg_list)


def init_arg(arg):
    # since visualize function is not strictly bound with other modules, only
    # check input here
    if arg.input_csv is None:
        log.error('Input is empty.')
        return None
    else:
        arg.input_csv = Path(arg.input_csv).absolute()
        if not arg.input_csv.exists():
            log.error('Input file does not exist.')
    if arg.og_type is None:
        log.error('Organelle type is empty.')
        return None
    if arg.ref_gb is None and arg.taxon is None:
        log.error('Reference is empty.')
        return None
    elif arg.ref_gb is not None:
        arg.ref_gb = Path(arg.ref_gb).absolute()
        if not arg.ref_gb.exists():
            log.error(f'{arg.ref_gb} does not exist.')
            return None
        else:
            log.info(f'Use {arg.ref_gb} as reference')
    elif arg.taxon is not None:
        arg.ref_gb = get_ref_gb(arg.taxon, arg.og_type)
    else:
        # if both arg.ref and arg.taxon is provided, use arg.ref
        pass
    arg.out = Path(arg.out).absolute()
    if arg.out.exists():
        log.warning(f'{arg.out} already exists. Overwriting.')
    # tobacco plastid genome
    arg.parts = {'LSC': arg.lsc, 'IRa': arg.ir, 'SSC': arg.ssc, 'IRb': arg.ir}
    arg.gb_size = (arg.parts['LSC'] + arg.parts['IRa'] + arg.parts['SSC'] +
                   arg.parts['IRb'])
    # total length - ir length
    arg.ir2_start = arg.lsc + arg.ir + arg.ssc
    return arg


def get_name(old):
    return '-'.join(old.split('-')[1:])


def get_name2(old):
    return old.split('-')[0]


def get_ref_gb(taxa: str, og_type: str) -> Path:
    log.info(f'Try to get reference of {taxa} {og_type} genome from Genbank')
    tmp_out = Path().cwd() / 'visualize_tmp'
    # Entrez.email = 'example.org'
    # if og_type == 'mt':
    #     query_str = 'mitochondrion[filter]'
    # else:
    #     query_str = '(plastid[filter] OR chloroplast[filter])'
    # query_str += ' AND Refseq[filter] AND "{taxa}"[organism]'
    # query_result = Entrez.read(Entrez.esearch(db='nuccore', term=query_str,
    #                                           usehistory='y'))
    # count = int(query_result['Count'])
    # if count == 0:
    #     log.error(f'{taxa} {og_type} genome not found in Genbank')
    #     return Path()
    # try:
    #     data = Entrez.efetch(db='nuccore',
    #                          webenv=query_result['WebEnv'],
    #                          query_key=query_result['QueryKey'],
    #                          rettype='acc',
    #                          retstart=0,
    #                          retmax=1)
    #     out_gb.write_text(data.read())
    # except Exception:
    #     log.critical(f'{taxa} {og_type} genome not found in Genbank')
    #     return Path()
    arg_str = (f' -taxon {taxa} -og {og_type} -out {tmp_out} -refseq yes '
               f'-count 1 -out_debug ')
    arg_result, _ = gb2fasta_main(arg_str)
    extend_gb = tmp_out / 'extend.gb'
    if arg_result.gb is not None and extend_gb.exists():
        log.info(f'Got {extend_gb} as reference genome for visualization')
    else:
        log.critical(f'Failed to get reference genome for visualization')
        return Path()
    return extend_gb


def percent_to_float(old: str) -> float:
    return float(old.rstrip('%')) / 100


# def read_csv(filename: Path) -> dict:
#     # Pandas is too heavy, use csv to implement pandas.read_csv
#     columns = defaultdict(list)
#     with open(filename, newline='') as f:
#         reader = csv.DictReader(f)
#         for row in reader:
#             for (k, v) in row.items():
#                 columns[k].append(v)
#     return columns


def preprocess_data(csv_file: Path, count_threshold) -> (pd.DataFrame, set):
    data_raw = pd.read_csv(csv_file)
    # spacer-atpB-rbcL to ('spacer', 'atpB-rbcL')
    data_raw['Name'] = data_raw.Loci.apply(get_name)
    data_raw['Loci'] = data_raw.Loci.apply(get_name2)
    data_raw['Gap_Ratio'] = data_raw['Gap_Ratio'].apply(percent_to_float)
    data_raw['Tree_Res'] = data_raw['Tree_Res'].apply(percent_to_float)
    data_raw['Total_GC'] = data_raw['Total_GC'].apply(percent_to_float)
    data_raw['Observed_Res'] = data_raw['Observed_Res'].apply(percent_to_float)
    data_raw2 = data_raw[data_raw.Samples >= count_threshold]
    # gene is duplicated with CDS/rna
    data_raw3 = data_raw2[data_raw2.Loci != 'gene']
    data_names_raw = data_raw3.Name.tolist()
    data_names = list()
    for i in data_names_raw:
        new = re.sub('trn(.|fM)_...', r'trn\1', i)
        data_names.append(new)
    names_set = set(data_names)
    return data_raw3, names_set


def draw_bar(track, x_, y_, w_, text, arg):
    max_y = max(y_)
    color = track_colors.pop()
    # ec='black', lw=0.5
    for x, y, w in zip(x_, y_, w_):
        # bar need x[list], y[list]
        # x should adjust by width
        track.bar([x - w / 2], [y], width=w * 0.95, color=color,
                  vmin=0, vmax=max_y)
    if arg.og_type == 'cp':
        t_pos = arg.gb_size - arg.parts['IRb'] / 2
        track.text(text, t_pos, size=8, color=color, adjust_rotation=True)


def visualize_main(arg_str=None):
    log.info('Running visualize module...')
    if arg_str is None:
        arg, other_args2 = parse_args()
    else:
        arg, other_args2 = parse_args(arg_str.split(' '))
    arg = init_arg(arg)
    if arg is None:
        log.info('Quit visualize module.')
        return None, other_args2

    # evaluation result table, remove abnormal genes
    long_label = 20
    # reference organelle genome genbank file, generated from OGU.gb2fasta
    gb = Genbank(arg.ref_gb)
    data_raw3, data_names_set = preprocess_data(arg.input_csv,
                                                arg.count_threshold)

    r = MyRadius()

    # genome gb file as template, no extra treatment
    if arg.og_type == 'cp':
        features = ('CDS', 'intron', 'tRNA', 'rRNA', 'spacer')
        feature_colors = '#FFD93D,#DDDDDD,#4D96FF,#9D96FF,#6BCB77'.split(',')
        circle_start = -260
        circle_end = circle_start + 360
    else:
        features = ('CDS', 'D-loop', 'tRNA', 'rRNA', 'spacer')
        feature_colors = '#FFD93D,#999999,#4D96FF,#9D96FF,#6BCB77'.split(',')
        circle_start = 0
        circle_end = 360
    circos = Circos(sectors={gb.name: gb.range_size},
                    sector2clockwise={gb.name: False},
                    start=circle_start, end=circle_end)
    sector = circos.sectors[0]

    r1 = r.get()
    feature_track = sector.add_track(r1)
    feature_track.axis(fc='#FFFFFF', ec='none')
    gene_labels = []
    gene_list = []
    last_gene = ''
    gene_set = set()
    ir_pos_list = set()
    # across_ir_ssc = {'ycf1', 'ndhF'}
    across_ir_ssc = set()
    for feat in gb.extract_features('gene'):
        start, end = int(str(feat.location.start)), int(str(feat.location.end))
        pos = (start + end) / 2
        label = feat.qualifiers.get('gene', ['??'])[0]
        if (label == '' or label.startswith('hypothetical')
                or label == last_gene):
            continue
        if len(label) > long_label:
            label = label[:long_label] + "..."
        last_gene = label
        gene_list.append(pos)
        gene_labels.append(label.rstrip("'"))
        if label in gene_set:
            ir_pos_list.add(pos)
        gene_set.add(label)
        # ycf1 or others span across ir
        if start < arg.ir2_start < end:
            across_ir_ssc.add(label)

    # print(across_ir_ssc)
    # Plot gene labels on outer position
    feature_track.xticks(
        gene_list,
        gene_labels,
        label_orientation='vertical',
        show_bottom_line=True,
        label_size=8,
        # line_kws=dict(ec="white"),
    )
    # extract names
    pos_list, labels = [], []
    last_name = ''
    for feature, color in zip(features, feature_colors):
        try:
            feature_track.genomic_features(
                gb.extract_features(feature),
                plotstyle='box',
                r_lim=r1,
                fc=color, alpha=1, lw=0.1, ec='black')
            for f in gb.extract_features(feature):
                start, end = int(str(f.location.start)), int(str(f.location.end))
                pos = end
                if feature == 'intron':
                    label = (f.qualifiers['gene'][0] + '.' +
                             f.qualifiers.get('number', '1')[0])
                elif feature != 'spacer':
                    label = f.qualifiers.get("gene", ["??"])[0]
                    if label == last_name:
                        continue
                else:
                    label = (f.qualifiers['upstream'][0] + '-'
                             + f.qualifiers['downstream'][0])
                last_name = label
                pos_list.append(pos)
                labels.append(label)
        except IndexError:
            log.info(f'Skip {feature}')

    labels_set = set(labels)
    intersection = data_names_set & labels_set
    clean_data = data_raw3.query('Name in @intersection')
    label_pos = list()
    label_pos_dict = dict()
    ir_second_label_pos = list()
    for key, value in zip(labels, pos_list):
        if key not in intersection:
            continue
        label_pos.append((key, value))
        # value_ = max(value, gene_end.get(key, 0))
        value_ = value
        if ((arg.og_type == 'cp' and value_ < arg.ir2_start
             or key in across_ir_ssc) or arg.og_type == 'mt'):
            label_pos_dict[key] = value
        else:
            ir_second_label_pos.append((key, value))

    # clean_data['Position'] = [label_pos_dict[x] for x in clean_data.Name]
    position_list = []
    for x in clean_data.Name:
        if x in label_pos_dict:
            position_list.append(label_pos_dict[x])
        else:
            position_list.append(0)
            log.info(f'Skip {x}')
    clean_data = clean_data.assign( Position=position_list)
    clean_data = clean_data.sort_values(by=['Position'])
    clean_pos = clean_data.Position.tolist()
    widths = list()
    widths.append(clean_pos[0])
    for i in range(1, len(clean_pos)):
        w = clean_pos[i] - clean_pos[i - 1]
        widths.append(w)

    for col, text in colname_text:
        track = sector.add_track(r.get(), r_pad_ratio=0.1)
        track.axis()
        data = clean_data[col].tolist()
        draw_bar(track, clean_pos, data, widths, text, arg)

    if arg.og_type == 'cp':
        # draw quadripartite structure for plastid
        part_track = sector.add_track(r.get(True))
        start = 0
        for part in arg.parts.keys():
            part_track.rect(start, start + arg.parts[part], fc=part_color.pop(),
                            lw=0.5)
            part_track.text(part, (start + start + arg.parts[part]) / 2,
                            color='white', size=8)
            start = start + arg.parts[part] + 1
        circos.link(
            (gb.name, arg.parts['LSC'], arg.parts['LSC'] + arg.parts['IRa']),
            (gb.name, arg.gb_size - arg.parts['IRb'], arg.gb_size),
            color='green', alpha=0.3)

    if arg.og_type == 'mt':
        fig = circos.plotfig(figsize=(14, 14))
    else:
        fig = circos.plotfig(figsize=(12, 12))

    # add legends
    rect_handles = []
    for name, color in zip(features, feature_colors):
        rect_handles.append(Patch(color=color, label=name))
    if arg.og_type == 'mt':
        # rect_handles.append(Patch(color='#999999', label='D-loop'))
        for name, color in zip(colname_text, track_colors2):
            rect_handles.append(Patch(color=color, label=name[1]))
        _ = circos.ax.legend(handles=rect_handles, bbox_to_anchor=(1, 0),
                             loc='lower center', fontsize=10, title='Types',
                             ncol=3)
    else:
        _ = circos.ax.legend(handles=rect_handles,
                             # bbox_to_anchor=(0.5, 0.5),
                             loc='lower right', fontsize=10, title='Types',
                             ncol=2)
    fig.savefig(arg.out)
    return arg, other_args2


if __name__ == '__main__':
    visualize_main()
