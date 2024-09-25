import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from pycirclize import Circos
from pycirclize.parser import Genbank

from OGU.global_vars import log

count_threshold = 100
# Evaluation.csv -> xlsx -> first column Type, second column name
# total length - ir length
irb_start = 131597
parts = {'LSC': 86684, 'IRa': 25339, 'SSC': 18482, 'IRb': 25339}
feature_colors = '#FFD93D,#DDDDDD,#4D96FF,#9D96FF,#6BCB77'.split(',')
part_color = '#6BCB77,#4D96FF,#6BCB77,#FF6B6B'.split(',')
# yellow,yellow/gray,blue,blue,green?
track_colors = list(reversed('#F7931B,#66B3FF,#338CFF,#3F51B5,#FB9883,#CC5500,'
                             '#FF6347,#36454F,#778877'.split(',')))
features = ('CDS', 'intron', 'tRNA', 'rRNA', 'spacer')
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
    arg.add_argument('-input', help='Evaluation result tabel, csv file')
    arg.add_argument('-ref', help='reference genome genbank file')
    arg.add_argument('-taxa', help='Taxa name, recommend use family name')
    arg.add_argument('-type', choices=('cp', 'mt'), dest='og_type',
                     help='Type of organelle (mt:mitochondria or cp:plastid')
    arg.add_argument('n_seqs', type=int, dest='count_threshold', default=10,
                     help='Minimum number of sequences per gene/fragment')
    # use tobacco as default
    arg.add_argument('-lsc', type=int, default=86684,
                     help='LSC size of plastid genome')
    arg.add_argument('-ssc', type=int, default=18482,
                     help='SSC size of plastid genome')
    arg.add_argument('-ir', type=int, default=25339,
                     help='IR size of plastid genome')
    arg.add_argument('-out', help='Output file', default='figure.pdf')
    if arg_list is None:
        return arg.parse_known_args()
    else:
        return arg.parse_known_args(arg_list)


def init_args(arg):
    # since visualize function is not strictly bound with other modules, only
    # check input here
    if arg.input is None:
        log.error('Input is empty.')
        return None
    else:
        arg.input = Path(arg.input).absolute()
        if not arg.input.exists():
            log.error('Input file does not exist.')
    if arg.og_type is None:
        log.error('Organelle type is empty.')
        return None
    if arg.ref is None and arg.taxa is None:
        log.error('Reference is empty.')
        return None
    elif arg.ref is not None:
        arg.ref = Path(arg.ref).absolute()
        if not arg.ref.exists():
            log.error(f'{arg.ref} does not exist.')
            return None
        else:
            log.info(f'Use {arg.ref} as reference')
    elif arg.taxa is not None:
        arg.ref = get_ref_gb(arg.taxa, arg.og_type)
    else:
        # if both arg.ref and arg.taxa is provided, use arg.ref
        pass
    arg.out = Path(arg.out).absolute()
    if arg.out.exists():
        log.warning(f'{arg.out} already exists. Overwriting.')
    return arg


def get_ref_gb(taxa: str, og_type: str) -> Path:
    log.info(f'Try to get reference of {taxa} {og_type} genome from Genbank')
    pass


def get_name(old):
    return '-'.join(old.split('-')[1:])


def get_name2(old):
    return old.split('-')[0]


def percent_to_float(old: str) -> float:
    return float(old.rstrip('%')) / 100


def read_csv(filename: Path) -> dict:
    # Pandas is too heavy, use csv to implement pandas.read_csv
    columns = defaultdict(list)
    with open(filename, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for (k, v) in row.items():
                columns[k].append(v)
    return columns


def x():
    #%%
    # read data
    # evaluation result table, remove abnormal genes
    filename = 'Evaluation.csv'
    # reference organelle genome genbank file, generated from OGU.gb2fasta
    gb_file = r'tobacco_extend.gb'
    gb = Genbank(gb_file)
    data_raw = pd.read_csv(filename)
    data_raw['Name'] = data_raw.Loci.apply(get_name)
    data_raw['Loci'] = data_raw.Loci.apply(get_name2)
    data_raw['Gap_Ratio'] = data_raw['Gap_Ratio'].apply(percent_to_float)
    data_raw['Tree_Res'] = data_raw['Tree_Res'].apply(percent_to_float)
    data_raw['Total_GC'] = data_raw['Total_GC'].apply(percent_to_float)
    data_raw['Observed_Res'] = data_raw['Observed_Res'].apply(percent_to_float)
    data_raw2 = data_raw[data_raw.Samples >= count_threshold]
    data = data_raw2[data_raw2.Loci != 'gene']
    data_names_raw = data.Name.tolist()
    data_names = list()
    for i in data_names_raw:
        new = re.sub('trn(.|fM)_...', r'trn\1', i)
        data_names.append(new)

    #%%
    r = MyRadius()

    #%%
    fig = plt.figure(figsize=(10, 10))
    # circos = Circos(sectors, space=0.5)
    # genome gb file as template, no extra treatmetn
    circos = Circos(sectors={gb.name: gb.range_size},
                    sector2clockwise={gb.name: False},
                    start=-260, end=-260 + 360)
    sector = circos.sectors[0]

    r1 = r.get()
    feature_track = sector.add_track(r1)
    feature_track.axis(fc="#FFFFFF", ec="none")
    gene_labels = []
    gene_list = []
    last_gene = ''
    gene_set = set()
    ir_pos_list = set()
    for feat in gb.extract_features("gene"):
        start, end = int(str(feat.location.start)), int(str(feat.location.end))
        pos = (start + end) / 2
        label = feat.qualifiers.get("gene", ["??"])[0]
        if label == "" or label.startswith("hypothetical") or label == last_gene:
            continue
        if len(label) > 20:
            label = label[:20] + "..."
        last_gene = label
        gene_list.append(pos)
        gene_labels.append(label)
        if label in gene_set:
            ir_pos_list.add(pos)
        gene_set.add(label)

    # Plot CDS product labels on outer position
    feature_track.xticks(
        gene_list,
        gene_labels,
        label_orientation="vertical",
        show_bottom_line=True,
        label_size=6,
        #line_kws=dict(ec="white"),
    )
    # extract names
    pos_list, labels = [], []
    last_name = ''
    for feature, color in zip(features, feature_colors):
        feature_track.genomic_features(
            gb.extract_features(feature),
            plotstyle="box",
            r_lim=r1,
            fc=color, alpha=1, lw=0.1, ec='black')
        for f in gb.extract_features(feature):
            start, end = int(str(f.location.start)), int(str(f.location.end))
            pos = end
            if feature == 'intron':
                label = f.qualifiers['gene'][0] + '.' + f.qualifiers['number'][0]
            elif feature != 'spacer':
                label = f.qualifiers.get("gene", ["??"])[0]
                if label == last_name:
                    continue
            else:
                label = f.qualifiers['upstream'][0] + '-' + f.qualifiers['downstream'][0]
            last_name = label
            pos_list.append(pos)
            labels.append(label)

    labels_set = set(labels)
    data_names_set = set(data_names)
    #a_b = sorted(list(data_names_set-labels_set))
    #print('data_names_set-labels_set', a_b, len(a_b))
    #print('\n'*3)
    #b_a = sorted(list(labels_set-data_names_set))
    #print('labels_set-data_names_set', b_a, len(b_a))
    intersection = data_names_set & labels_set
    clean_data = data.query('Name in @intersection')
    label_pos = list()
    label_pos_dict = dict()
    ir_second_label_pos = list()
    for key, value in zip(labels, pos_list):
        if key not in intersection:
            continue
        label_pos.append((key, value))
        # ycf1, irb start
        if value < irb_start:
            label_pos_dict[key] = value
        else:
            ir_second_label_pos.append((key, value))

    clean_data['Position'] = [label_pos_dict[x] for x in clean_data.Name]
    clean_data = clean_data.sort_values(by=['Position'])
    clean_pos = clean_data.Position.tolist()
    widths = list()
    widths.append(clean_pos[0])
    for i in range(1, len(clean_pos)):
        w = clean_pos[i] - clean_pos[i - 1]
        widths.append(w)

    t_pos = gb.range_size - parts['IRb'] / 2


    def draw_bar(track, x_, y_, w_, text):
        max_y = max(y_)
        color = track_colors.pop()
        for x, y, w in zip(x_, y_, w_):
            # bar need x[list], y[list]
            # x should adjust by width
            track.bar([x - w / 2], [y], width=w * 0.95, color=color,
                      vmin=0, vmax=max_y)
        track.text(text, t_pos, size=8, color=color, adjust_rotation=True)

        # ec='black', lw=0.5


    for col, text in colname_text:
        track = sector.add_track(r.get(), r_pad_ratio=0.1)
        track.axis()
        data = clean_data[col].tolist()
        draw_bar(track, clean_pos, data, widths, text)

    part_track = sector.add_track(r.get(True))
    start = 0
    for part in parts.keys():
        part_track.rect(start, start + parts[part], fc=part_color.pop(), lw=0.5)
        part_track.text(part, (start + start + parts[part]) / 2,
                        color='white', size=8)
        start = start + parts[part] + 1
    circos.link((gb.name, parts['LSC'], parts['LSC'] + parts['IRa']),
                (gb.name, gb.range_size - parts['IRb'], gb.range_size),
                color='green', alpha=0.3)

    # print(clean_data.columns)

    fig = circos.plotfig()

    rect_handles = []
    rect_colors = feature_colors
    for name, color in zip(features, rect_colors):
        rect_handles.append(Patch(color=color, label=name))
    _ = circos.ax.legend(handles=rect_handles, # bbox_to_anchor=(0.5, 0.5),
                         loc='lower right', fontsize=8, title="Types", ncol=2)
    #%%
    fig.savefig('out.pdf')


def visualize_main():
    pass


if __name__ == '__main__':
    visualize_main()
