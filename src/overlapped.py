"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
import numpy as np
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import Span, BoxAnnotation, HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, CustomJS
from bokeh.colors import RGB
import pyslow5
import copy
import argparse
import re
import logging

BASE_SHIFT = 5
BASE_LIMIT = 1000
FRONT_CLIP = 50
SIG_PLOT_LENGTH = 8000
TRIM_OFFSET = 0
base_color_map = {'A': 'limegreen', 'C': 'blue', 'T': 'red', 'G': 'orange'}

# # importing the modules
# from bokeh.plotting import figure, output_file, show
# # file to save the model
# output_file("gfg.html")
# # instantiating the figure object
# graph = figure(title="Bokeh Multi Line Graph")
# # the points to be plotted
# xs = [[1, 2, 3, 4, 5], [-4, -2, 0, 2, 4]]
# ys = [[5, 3, 8, 0], [5, -4, 10, -2, 5]]
# # plotting the graph
# graph.multi_line(xs, ys)
# # displaying the model
# show(graph)
# exit(0)

parser = argparse.ArgumentParser()
   
parser.add_argument('-f', '--fastq', required=False, help="fastq file")
parser.add_argument('-r', '--read_id', required=False, help="read id")
parser.add_argument('--trim_offset', required=False, help="signal trim offset")
parser.add_argument('--plot_signal_limit', required=False, help="signal plot length")
parser.add_argument('-s', '--slow5', required=True, help="slow5 file")
parser.add_argument('--second_slow5', required=True, help="second slow5 file")
parser.add_argument('--second_fastq', required=False, help="second fastq file")
parser.add_argument('-a', '--alignment', required=False, help="alignment file")
parser.add_argument('-b', '--base_shift', required=False, default=BASE_SHIFT, help="base shift")
parser.add_argument('-o', '--output', required=True, help="output file (html)")

args = parser.parse_args()

flag_read_id = 0
if args.read_id:
    flag_read_id = 1
    print(f'read_id: {args.read_id}')
    read_id = args.read_id

trim_offset = TRIM_OFFSET
if args.trim_offset:
    trim_offset = args.trim_offset

plot_signal_limit = SIG_PLOT_LENGTH
if args.plot_signal_limit:
    plot_signal_limit = args.plot_signal_limit

flag_fastq = 0
if args.fastq:
    flag_fastq = 1
    print(f'fastq file: {args.fastq}')
print(f'signal file: {args.slow5}')
flag_alignment = 0
if args.alignment:
    flag_alignment = 1
    print(f'alignment file: {args.alignment}')
print(f'base_shift: {args.base_shift}')
print(f'output file: {args.output}')

if flag_alignment == 1 and flag_fastq == 0:
    logging.error("please provide both fastq and alignment files")
    exit()

if flag_fastq == 1:
    # read fastq file
    fastq_file = open(args.fastq, 'r')
    read_id = fastq_file.readline().split()[0][1:]
    fastq_seq = fastq_file.readline().rstrip()
    if BASE_LIMIT > len(fastq_seq):
        BASE_LIMIT = len(fastq_seq)
    fastq_file.close()

if flag_alignment == 1:
    # read alignment file
    align_file = open(args.alignment, 'r')
    alignment = align_file.readline().split()
    trim_offset = int(alignment[2])
    moves_string = alignment[14]
    align_file.close()

# open signal file
x = []
x_real = []
y = [0]

s5 = pyslow5.Open(args.slow5, 'r')
read = s5.get_read(read_id, pA=True, aux=["read_number", "start_mux"])
if read is not None:
    print("read_id:", read['read_id'])
    print("len_raw_signal:", read['len_raw_signal'])
    # x = list(range(1,read['len_raw_signal']+1))
    # y = read['signal']

    chunk_offset = 0
    if trim_offset == 0:
        FRONT_CLIP = 0
    start_index = trim_offset + chunk_offset - FRONT_CLIP
    if read['len_raw_signal'] < start_index + SIG_PLOT_LENGTH:
        SIG_PLOT_LENGTH = read['len_raw_signal'] - start_index
        end_index = read['len_raw_signal']
    else:
        end_index = start_index + SIG_PLOT_LENGTH
    shift = end_index

    x = list(range(0, end_index-start_index + 1))
    x_real = list(range(start_index, end_index + 1))
    y.extend(read['signal'][start_index:end_index])
s5.close()
# set output to static HTML file
output_file(filename=args.output, title=read_id)

plot_title = f'{read_id}:{start_index}-{end_index}-{trim_offset}-{chunk_offset}'
tools_to_show = 'hover,box_zoom,pan,save,wheel_zoom'
p = figure(title=plot_title,
		   x_axis_label='signal index',
		   y_axis_label='signal value',
		   sizing_mode="stretch_width",
		   height=300,
		   output_backend="webgl",
		   x_range=(0, 750),
		   tools=tools_to_show)
    # tooltips=tool_tips)

p.toolbar.active_scroll = p.select_one(WheelZoomTool)

if flag_alignment == 1:
    base_x = []
    base_y = []
    base_label = []
    base_count = 0
    location_plot = trim_offset+chunk_offset - start_index
    previous_location = location_plot
    # draw moves
    moves_string = re.sub('ss:Z:', '', moves_string)
    moves_string = re.sub('D', 'D,', moves_string)
    moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
    # print(moves_string)
    moves = re.split(r',+', moves_string)

    vlines = []
    BASE_SHIFT = int(args.base_shift)
    base_count = BASE_SHIFT

    for i in moves:
        previous_location = location_plot
        n_samples = 0
        if 'D' in i:
            i = re.sub('D', '', i)
            n_samples = int(i)
            prev_loc = previous_location
            for j in range(0, n_samples):
                base = fastq_seq[base_count]
                base_box = BoxAnnotation(left=prev_loc, right=prev_loc+5, fill_alpha=0.2, fill_color='white')
                p.add_layout(base_box)

                base_x.append(prev_loc)
                base_y.append(115)
                label = str(base) + "\t" + str(base_count + 1)
                base_label.append(label)

                prev_loc = prev_loc + 5
                x_end = x[-1]
                x = x + list(range(x_end, x_end+5))
                base_count = base_count + 1
            location_plot = prev_loc
            z = np.concatenate((y[:previous_location], [0] * n_samples * 5), axis=0)
            y = np.concatenate((z, y[previous_location:]), axis=0)
            # y = y[:previous_location] + y[previous_location:]

        elif 'I' in i:
            i = re.sub('I', '', i)
            n_samples = int(i)
            location_plot = location_plot + n_samples

            vline = Span(location=location_plot, dimension='height', line_color='red', line_width=1)
            vlines.append(vline)

        else:
            n_samples = int(i)
            location_plot = location_plot + n_samples

            base = fastq_seq[base_count]
            base_box = BoxAnnotation(left=previous_location, right=location_plot, fill_alpha=0.2, fill_color=base_color_map[base])
            p.add_layout(base_box)

            vline = Span(location=location_plot, dimension='height', line_color='red', line_width=1)
            vlines.append(vline)

            base_x.append(previous_location)
            base_y.append(115)
            label = str(base) + "\t" + str(base_count + 1)
            base_label.append(label)
            base_count = base_count + 1

        if base_count == BASE_LIMIT:
            break
        if location_plot > SIG_PLOT_LENGTH:
            break

    p.renderers.extend(vlines)

    base_annotation = ColumnDataSource(data=dict(base_x=base_x,
                                                 base_y=base_y,
                                                 base_label=base_label))

    base_annotation_labels = LabelSet(x='base_x', y='base_y', text='base_label',
                                      x_offset=5, y_offset=5, source=base_annotation, render_mode='canvas',
                                      text_font_size="7pt")

    p.add_layout(base_annotation_labels)

    plot_signal_limit = location_plot + 10

y_second = [0]


if flag_fastq == 1:
    # read fastq file
    fastq_file = open(args.second_fastq, 'r')
    read_id = fastq_file.readline().split()[0][1:]
    fastq_seq = fastq_file.readline().rstrip()
    if BASE_LIMIT > len(fastq_seq):
        BASE_LIMIT = len(fastq_seq)
    fastq_file.close()

s5 = pyslow5.Open(args.second_slow5, 'r')
read = s5.get_read(read_id, pA=True)
if read is not None:
    print("read_id:", read['read_id'])
    print("len_raw_signal:", read['len_raw_signal'])
    # x = list(range(1,read['len_raw_signal']+1))
    # y = read['signal']
    chunk_offset = 0
    if trim_offset == 0:
        FRONT_CLIP = 0
    start_index = trim_offset + chunk_offset - FRONT_CLIP
    if read['len_raw_signal'] < start_index + SIG_PLOT_LENGTH:
        SIG_PLOT_LENGTH = read['len_raw_signal'] - start_index
        end_index = read['len_raw_signal']
    else:
        end_index = start_index + SIG_PLOT_LENGTH
    shift = end_index

    y_second.extend(read['signal'][start_index:end_index])
s5.close()

SIGNAL_SHIFT = 15
y_second = [value + SIGNAL_SHIFT for value in y_second]

p.line(x, y, line_width=2, color='black', legend_label='nanopolsih')
p.line(x, y_second, line_width=1, color='grey', legend_label='a14')
# add a circle renderer with a size, color, and alpha
# p.circle(x[:plot_signal_limit], y[:plot_signal_limit], size=2, color="red", alpha=0.5)

# show the tooltip
hover = p.select(dict(type=HoverTool))
hover.tooltips = [("x", "@x_real"), ("y", "$y")]
hover.mode = 'mouse'

save(p)
