import numpy as np
from bokeh.io import output_file, show
from bokeh.layouts import gridplot
from bokeh.plotting import figure, save
from bokeh.models import Span, BoxAnnotation, HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, CustomJS
import argparse

FILLING_SIZE = 30

parser = argparse.ArgumentParser()

parser.add_argument('-k', '--kmer', required=True, help="kmer")
parser.add_argument('-d', '--kmer_dump', required=True, help="kmer dump file")
parser.add_argument('-o', '--output', required=True, help="output file (html)")

args = parser.parse_args()
print(f'kmer: {args.kmer}')
print(f'kmer_dump: {args.kmer_dump}')
print(f'output file: {args.output}')
output_file(filename=args.output, title=args.kmer)

# read kmer_dump file
kmer_dump_file = open(args.kmer_dump, 'r')
kmer_dump = kmer_dump_file.readline().split(';')[:-1]
kmer_dump_file.close()

n_kmer_dump = len(kmer_dump)
print(f'number of plots : {n_kmer_dump}')

filling = [0] * FILLING_SIZE

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

x = []
y = []

start_count = 0

for i in kmer_dump:

    y_sig = i.split(',')
    len_i = len(y_sig)
    x_sig = list(range(start_count,  start_count + len_i + FILLING_SIZE))
    start_count = start_count + len_i + FILLING_SIZE

    x = x + x_sig
    y = y + y_sig + filling

# print(x)
# print(y)
# print(len(x))
# print(len(y))

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
