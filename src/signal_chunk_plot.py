import numpy as np
from bokeh.io import output_file, show
from bokeh.layouts import gridplot
from bokeh.plotting import figure, save
from bokeh.models import Span, BoxAnnotation, HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, CustomJS
import argparse

from numpy import array, linspace
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema

FILLING_SIZE = 30

def main():

    flag_fastq_on = 0

    parser = argparse.ArgumentParser()

    parser.add_argument('-k', '--kmer', required=True, help="kmer")
    parser.add_argument('-d', '--kmer_dump', required=True, help="kmer dump file")
    parser.add_argument('-m', '--signal_print_margin', required=True, help="signal print margin used in poregen")
    parser.add_argument('-o', '--output', required=True, help="output file (html)")
    parser.add_argument('-r', '--fastq', required=False, help="fastq file (optional)")
    parser.add_argument('-n', '--ngb', required=False, help="no. of bases to the right (optional)")

    args = parser.parse_args()
    print(f'kmer: {args.kmer}')
    print(f'kmer_dump: {args.kmer_dump}')
    print(f'signal_print_margin: {args.signal_print_margin}')
    print(f'output file: {args.output}')

    kmer_size = len(args.kmer)
    kmer_neigbours_right = []

    num_neighbours = 1
    if args.ngb:
        num_neighbours = int(args.ngb)

    if args.fastq:
        print(f'fastq file: {args.fastq}')
        flag_fastq_on = 1
        fastq_file = open(args.fastq, 'r')
        fastq_file.readline()
        fastq_seq = fastq_file.readline().rstrip()
        shift_count = 0
        fastq_len = len(fastq_seq)
        for s in fastq_seq[kmer_size:]:
            kmer_ = fastq_seq[shift_count:shift_count+kmer_size]
            if kmer_ == args.kmer and shift_count+kmer_size+num_neighbours < fastq_len:
                bases_to_right = fastq_seq[shift_count+kmer_size+1:shift_count+kmer_size+1+num_neighbours]
                kmer_neigbours_right.append(bases_to_right)
            shift_count = shift_count + 1


    output_file(filename=args.output, title=args.kmer)

    # read kmer_dump file
    kmer_dump_file = open(args.kmer_dump, 'r')
    kmer_dump = kmer_dump_file.readline().replace(':', '').split(';')[:-1]
    kmer_dump_file.close()

    n_kmer_dump = len(kmer_dump)
    print(f'number of plots : {n_kmer_dump}')

    filling = [0] * FILLING_SIZE
    signal_print_margin = int(args.signal_print_margin)

    kmer_size=len(args.kmer)

    plot_title = f'{kmer_size}-mer:{args.kmer}\t\t\tno. of plots:{n_kmer_dump}\t\t\tgap between siglets:{FILLING_SIZE}'
    tools_to_show = 'hover,box_zoom,pan,save,wheel_zoom'
    p = figure(title=plot_title,
               y_axis_label='signal value',
               sizing_mode="stretch_width",
               height=300,
               output_backend="webgl",
               x_range=(0, 750),
               tools=tools_to_show)
        # tooltips=tool_tips)

    p.toolbar.active_scroll = p.select_one(WheelZoomTool)

    x = [] # x values: margin, annotation, gap
    y = [] # y values: margin, annoration, gap

    neighbour_x = []
    neighbour_y = []

    start_count = 0

    guppy_y_all = [] # list to contain all annotated valued

    for i in kmer_dump:

        # y_sig = i.split(',')
        y_sig = [float(i) for i in i.split(',')]
        len_i = len(y_sig)
        x_sig = list(range(start_count,  start_count + len_i + FILLING_SIZE))
        start_count = start_count + len_i + FILLING_SIZE

        x = x + x_sig
        y = y + y_sig + filling

        guppy_x = x_sig[signal_print_margin:-(signal_print_margin+FILLING_SIZE)]
        guppy_y = y_sig[signal_print_margin:-signal_print_margin]
        guppy_y_all = guppy_y_all + guppy_y
        p.circle(guppy_x, guppy_y, size=5, color="red", alpha=0.5)
        calculate_current_mean = method_0(guppy_y)
        current_line = [calculate_current_mean for j in guppy_x]
        p.line(guppy_x, current_line, color='black', line_width=2, line_dash='dotdash')

        if flag_fastq_on:
            x_value = guppy_x[-1]
            neighbour_x.append(x_value)
            neighbour_y.append(120)

    # print(guppy_y_all)
    cluster_median = method_1(guppy_y_all)
    for median in cluster_median:
        p.line(x, median, color='grey', line_width=2, line_dash='dotdash')

    source_1 = ColumnDataSource(data=dict(x=neighbour_x,
                                        y=neighbour_y,
                                        bases=kmer_neigbours_right))
    labels = LabelSet(x='x', y='y', text='bases',
                      x_offset=5, y_offset=5, source=source_1)
    p.add_layout(labels)
    # print(x)
    # print(y)
    # print(len(x))
    # print(len(y))
    # p.circle(x, y, size=5, color="red", alpha=0.5)

    source = ColumnDataSource(data=dict(
        x=x,
        y=y,
    ))
    p.line('x', 'y', line_width=2, source=source)

    # show the tooltip
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [("x", "$x"), ("y", "$y")]
    hover.mode = 'mouse'

    save(p)


def method_0(y):
    med = np.median(y)
    return med

def method_1(y):
    cluster_median = []
    a = np.array(y).reshape(-1, 1)
    # a = samples.reshape(-1,1)
    kde = KernelDensity(kernel='gaussian', bandwidth=3).fit(a)
    s = linspace(0, 100)
    e = kde.score_samples(s.reshape(-1, 1))
    mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]

    if len(mi) == 0:
        print("one cluster")
        cluster_median.append(np.median(y))
    elif len(mi) == 1:
        print("two clusters")
        cluster = a[a < s[mi][0]]
        cluster_median.append(np.median(cluster))
        cluster = a[a >= s[mi][-1]]
        cluster_median.append(np.median(cluster))
    else:
        print("more than two clusters")
        cluster_median.append(np.median(y))

    return cluster_median

main()