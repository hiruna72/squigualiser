"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import BoxAnnotation, HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, Segment, Toggle, Range1d, FreehandDrawTool
import re

PLOT_X_RANGE = 300
PLOT_HEIGHT = 600

def adjust_before_plotting(ref_seq_len, signal_tuple, region_tuple, sig_algn_data, fasta_seq):
    if sig_algn_data["data_is_rna"]:
        ref_region_start_diff = region_tuple[1] - region_tuple[3]
    else:
        ref_region_start_diff = region_tuple[2] + 1 - region_tuple[0]
    ref_region_end_diff = region_tuple[2] + ref_seq_len - region_tuple[1]
    # print("ref_region_start_diff: " + str(ref_region_start_diff))
    moves = sig_algn_data['ss']
    if ref_region_start_diff < 0:
        # we have to remove part of the sequence until we get to the ref start
        count_bases = 0
        eat_signal = 0
        count_moves = 0
        prev_move = None
        updated_move = []
        for i in moves:
            if count_bases == abs(ref_region_start_diff):
                break
            if count_bases > abs(ref_region_start_diff):
                if not prev_move.find('D'):
                    raise Exception("Error: a deletion move was expected. incorrect implementation. Please report with a minimal reproducible test")
                updated_move = ["{}D".format(count_bases-abs(ref_region_start_diff))]
                break
            count_moves += 1

            if 'D' in i:
                i = re.sub('D', '', i)
                count_bases += int(i)
                # print(i+" D "+str(int(sam_record.reference_start) + 1 + count_bases))
            elif 'I' in i:
                i = re.sub('I', '', i)
                eat_signal += int(i)
                # print(i+" I "+str(int(sam_record.reference_start) + 1 + count_bases))
            else:
                eat_signal += int(i)
                count_bases += 1
            prev_move = i
        moves = updated_move + moves[count_moves:]
        x = signal_tuple[0][:-eat_signal]
        x_real = signal_tuple[1][eat_signal:]
        y = signal_tuple[2][eat_signal:]
        signal_tuple = (x, x_real, y)
        sig_algn_data['ss'] = moves
    return signal_tuple, region_tuple, sig_algn_data, fasta_seq

def create_figure(args):
    y_axis_label = "signal value (raw)"
    if args.no_pa:
        y_axis_label = "signal value (pA)"
    tools_to_show = 'hover,box_zoom,pan,reset,save,wheel_zoom,zoom_in,zoom_out'
    p_default = figure(x_axis_label='signal index',
               y_axis_label=y_axis_label,
               sizing_mode="stretch_width",
               height=PLOT_HEIGHT,
               output_backend="webgl",
               x_range=(0, PLOT_X_RANGE),
               tools=tools_to_show,
               toolbar_location="below")
    # tooltips=tool_tips)

    p_default.select(dict(type=WheelZoomTool)).maintain_focus = False
    p_default.toolbar.active_scroll = p_default.select_one(WheelZoomTool)
    p_default.toolbar.logo = None
    return p_default

