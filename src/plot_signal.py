"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
import bokeh
import numpy as np
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import BoxAnnotation, HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, Segment, Toggle, Range1d, FreehandDrawTool, CustomJS, FixedTicker, Spacer, Div
from bokeh.layouts import row, column
from bokeh.colors import RGB
from bokeh.io import export_svg, export_svgs
import pyslow5
import argparse
import re
import os
from src import plot_utils

DEFAULT_KMER_SIZE = 9
BASE_LIMIT = 1000
SIG_PLOT_LENGTH = 20000
DEFAULT_STRIDE = 5
PLOT_X_RANGE = 300
PLOT_HEIGHT = 600
PLOT_BASE_SHIFT = 0
PLOT_X_PADDING = 100
PLOT_LIMIT = 1000

DEFAULT_NUM_BED_COLS = 3
DEFAULT_BED_ANNOTATION_COLOR = (75, 126, 246)
BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK = range(10)
READ_ID, LEN_RAW_SIGNAL, START_RAW, END_RAW, STRAND, SEQUENCE_ID, LEN_KMER, START_KMER, END_KMER, MATCHES, LEN_KMER, MAPQ = range(12)
SI_START_RAW, SI_END_RAW, SI_START_KMER, SI_END_KMER = range(4)
BED_CHROM, BED_CHROM_START, BED_CHROM_END, BED_NAME, BED_SCORE, BED_STRAND, BED_THICK_START, BED_THICK_END, BED_ITEM_RGB, BED_BLOCK_COUNT, BED_BLOCK_SIZES, BLOCK_STARTS = range(12)

def plot_function(args, p, read_id, signal_tuple, draw_data):
    x = signal_tuple[0]
    x_real = signal_tuple[1]
    y = signal_tuple[2]

    y_min = draw_data['y_min']
    y_max = draw_data['y_max']

    source = ColumnDataSource(data=dict(x=x, y=y, x_real=x_real))

    p.line('x', 'y', name="sig_plot_line", line_width=2, source=source)
    # add a circle renderer with a size, color, and alpha
    sample_labels = p.circle(x, y, radius=draw_data["point_size"], color="red", alpha=0.5, visible=draw_data['no_samples'])
    toggle_samples = Toggle(label="samples", button_type="default", active=draw_data['no_samples'], height=30, width=60)
    toggle_samples.js_link('active', sample_labels, 'visible')

    # show the tooltip
    hover = p.select(dict(type=HoverTool))
    hover.renderers = p.select(name="sig_plot_line")
    hover.tooltips = [("x", "@x_real"), ("y", "@y")]
    hover.mode = 'mouse'

    p.title = args.tag_name

    location_plot = len(x_real)
    if location_plot > (y_max - y_min):
        if location_plot > draw_data['xrange']:
            p.x_range = Range1d(0, draw_data['xrange'], bounds=(-1*PLOT_X_PADDING, location_plot+PLOT_X_PADDING))
        else:
            p.x_range = Range1d(0, location_plot, bounds=(-1*PLOT_X_PADDING, location_plot+PLOT_X_PADDING))

    renderer = p.multi_line([[1, 1]], [[1, 1]], line_width=4, alpha=0.4, color='black')
    draw_tool = FreehandDrawTool(renderers=[renderer], num_objects=50)
    p.add_tools(draw_tool)

    
    message_browser = Div(text=f"Bokeh version: {bokeh.__version__} (Google Chrome is recommended)", width=400, height=30)

    layout_ = p, row(toggle_samples, Spacer(width=10)), column(message_browser)
    return layout_

def run(args):
    
    if args.read_id != "":
        args.plot_limit = 1


    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    
    if args.save_svg:
        plot_utils.svg_export_works(args.output_dir)  # will raise if it fails

    # open signal file
    print(f'signal file: {args.slow5}')
    s5 = pyslow5.Open(args.slow5, 'r')

    scaling_str = "no scaling"
    if args.sig_scale == "medmad" or args.sig_scale == "znorm":
        scaling_str = args.sig_scale
    elif not args.sig_scale == "":
        raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))

    scale_params = {}

    indt = "\t\t\t\t\t\t\t\t"
    draw_data = {}
    draw_data["point_size"] = args.point_size
    draw_data["no_samples"] = args.show_samples
    draw_data["xrange"] = args.xrange
    draw_data["sig_dir"] = "->"
    if args.reverse_signal:
        draw_data["sig_dir"] = "<-"

    start_index, end_index = args.region.split('-') if args.region else (1, args.sig_plot_limit)
    start_index = int(start_index) - 1  # convert to 0-based index
    end_index = int(end_index)  # keep it 1-based for slicing
    x, x_real, y = plot_utils.load_signal(args, args.read_id, s5, start_index, end_index, scale_params)

    # resolve -1 sentinel to actual signal length
    if end_index == -1:
        end_index = x_real[-1] if len(x_real) > 0 else start_index + 1
    if draw_data["xrange"] == -1:
        draw_data["xrange"] = len(x)

    if args.reverse_signal:
        # reverse the signal
        y = np.flip(y)
        # x = np.flip(x)
        x_real = np.flip(x_real)
        temp = start_index
        start_index = end_index - 1
        end_index = temp + 1

    signal_tuple = (x, x_real, y)
    draw_data['y_min'] = np.nanmin(signal_tuple[2])
    draw_data['y_max'] = np.nanmax(signal_tuple[2])

    output_file_name = args.output_dir + "/" + args.read_id + "_" + args.tag_name
    args.tag_name += f"{args.tag_name}{indt}{args.read_id}{indt}scale:{scaling_str}{indt}{draw_data['sig_dir']}{indt}region: [{start_index + 1}-{end_index}]" 
    
    args.fixed_width = False
    p = plot_utils.create_figure(args, plot_mode=0)

    layout_ = plot_function(args, p=p, read_id=args.read_id, signal_tuple=signal_tuple, draw_data=draw_data)

    if args.save_svg:
        output_file_name += ".svg"
        layout_[0].output_backend = "svg"
        export_svgs(layout_, filename=output_file_name)
    else:
        output_file_name += ".html"
        output_file(output_file_name, title=args.read_id)
        save(layout_)
    print(f'output file: {os.path.abspath(output_file_name)}')

    print(f'Bokeh version: {bokeh.__version__} (Google Chrome is recommended)')

    s5.close()

def argparser():
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    parser.add_argument('-r', '--read_id', required=True, type=str, default="", help="plot the read with read_id")
    parser.add_argument('-s', '--slow5', required=True, type=str, default="", help="slow5 file")
    parser.add_argument('--region', required=False, type=str, default="", help="[start-end] 1-based closed interval region to plot. eg:100-200.")
    parser.add_argument('--tag_name', required=False, type=str, default="", help="a tag name to easily identify the plot")
    parser.add_argument('--sig_scale', required=False, type=str, default="", help="plot the scaled signal. Supported scalings: [medmad, znorm, scaledpA]")
    parser.add_argument('--reverse_signal', required=False, action='store_true', help="reverse the signal")
    parser.add_argument('--no_pa', required=False, action='store_false', help="skip converting the signal to pA values")
    parser.add_argument('--remove_signal_outliers', required=False, action='store_true', help="remove signal outliers that are outside the raw value range [0, 2000]")
    parser.add_argument('--point_size', required=False, type=int, default=0.5, help="signal point radius [0.5]")
    parser.add_argument('--sig_plot_limit', required=False, type=int, default=-1, help="maximum number of signal samples to plot (-1 for full signal)")
    parser.add_argument('--show_samples', required=False, action='store_true', help="show sample points (default: hidden)")
    parser.add_argument('--save_svg', required=False, action='store_true', help="save as svg. tweak --region and --xrange to capture the necessary part of the plot")
    parser.add_argument('--xrange', required=False, type=int, default=-1, help="initial x range to display (-1 for full signal)")
    parser.add_argument('-o', '--output_dir', required=True, type=str, default="", help="output dir")
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)


