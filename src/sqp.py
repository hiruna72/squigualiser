"""
Signal to seQuence alignment Plot - sqp
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
from readpaf import parse_paf
from sys import stdin
from pyfaidx import Fasta
import os



base_color_map = {'A': 'limegreen', 'C': 'blue', 'T': 'red', 'G': 'orange'}

parser = argparse.ArgumentParser()
   
parser.add_argument('-f', '--fasta', required=True, help="fasta file")
parser.add_argument('-r', '--read_id', required=False, help="read id")
parser.add_argument('--plot_signal_limit', required=False, type=int, help="signal plot length")
parser.add_argument('-s', '--slow5', required=True, help="slow5 file")
parser.add_argument('-a', '--alignment', required=True, help="alignment file")
parser.add_argument('--point_size', required=False, type=int, default=5, help="signal point size [5]")
parser.add_argument('-o', '--output_dir', required=True, help="output dir")
def plot_function(read_id, output_file_name, signal_tuple, paf_record, fasta_sequence, base_limit, sig_plot_length):
    x = signal_tuple[0]
    x_real = signal_tuple[1]
    y = signal_tuple[2]

    plot_title = f'{read_id}  [signal_start_index,signal_end_index,signal_alignment_start_index:{start_index},{end_index},{paf_record.query_start}]  [seq_length,kmer_start_index,kmer_end_index:{paf_record.target_length},{paf_record.target_start},{paf_record.target_end}]'
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

    base_x = []
    base_y = []
    base_label = []
    base_count = 0
    location_plot = int(paf_record.query_start)
    previous_location = location_plot
    # draw moves
    moves_string = paf_record.tags['ss'][2]
    # moves_string = re.sub('ss:Z:', '', moves_string)
    moves_string = re.sub('D', 'D,', moves_string)
    moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
    # print(moves_string)
    moves = re.split(r',+', moves_string)

    vlines = []
    base_count = int(paf_record.target_start)

    for i in moves:
        previous_location = location_plot
        n_samples = 0
        if 'D' in i:
            i = re.sub('D', '', i)
            n_samples = int(i)
            prev_loc = previous_location
            for j in range(0, n_samples):
                base = fasta_sequence[base_count]
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

            base = fasta_sequence[base_count]
            base_box = BoxAnnotation(left=previous_location, right=location_plot, fill_alpha=0.2, fill_color=base_color_map[base])
            p.add_layout(base_box)

            vline = Span(location=location_plot, dimension='height', line_color='red', line_width=1)
            vlines.append(vline)

            base_x.append(previous_location)
            base_y.append(115)
            label = str(base) + "\t" + str(base_count + 1)
            base_label.append(label)
            base_count = base_count + 1

        if base_count == base_limit:
            break
        if location_plot > sig_plot_length:
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

    source = ColumnDataSource(data=dict(
        x=x[:plot_signal_limit],
        y=y[:plot_signal_limit],
        x_real=x_real[:plot_signal_limit],
    ))
    p.line('x', 'y', line_width=2, source=source)
    # add a circle renderer with a size, color, and alpha
    p.circle(x[:plot_signal_limit], y[:plot_signal_limit], size=args.point_size, color="red", alpha=0.5)

    # show the tooltip
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [("x", "@x_real"), ("y", "$y")]
    hover.mode = 'mouse'

    output_file(output_file_name, title=read_id)
    save(p)

args = parser.parse_args()

if args.fasta:
    print(f'fasta file: {args.fasta}')
print(f'signal file: {args.slow5}')
if args.alignment:
    print(f'alignment file: {args.alignment}')

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

# open signal file
s5 = pyslow5.Open(args.slow5, 'r')

with open(args.alignment, "r") as handle:
    fasta_reads = Fasta(args.fasta)
    for paf_record in parse_paf(handle):
        read_id = paf_record.query_name
        fasta_seq = fasta_reads.get_seq(name=read_id, start=1, end=int(paf_record.target_length)).seq
        output_file_name = args.output_dir+"/"+read_id+".html"
        print(f'output file: {output_file_name}')

        x = []
        x_real = []
        y = [0]
        base_limit = 1000
        sig_plot_length = 8000
        if args.plot_signal_limit:
            sig_plot_length = args.plot_signal_limit
        read = s5.get_read(read_id, pA=True, aux=["read_number", "start_mux"])
        if read is not None:
            print("read_id:", read['read_id'])
            print("len_raw_signal:", read['len_raw_signal'])
            start_index = 0
            if read['len_raw_signal'] < start_index + sig_plot_length:
                sig_plot_length = read['len_raw_signal'] - start_index
                end_index = read['len_raw_signal']
            else:
                end_index = start_index + sig_plot_length
            shift = end_index

            x = list(range(0, end_index - start_index + 1))
            x_real = list(range(start_index, end_index + 1))
            y.extend(read['signal'][start_index:end_index])
        signal_tuple = (x, x_real, y)

        plot_function(read_id=read_id, output_file_name=output_file_name, signal_tuple=signal_tuple, paf_record=paf_record, fasta_sequence=fasta_seq, base_limit=base_limit, sig_plot_length=sig_plot_length)

s5.close()
