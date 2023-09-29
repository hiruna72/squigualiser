"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import BoxAnnotation, HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, Segment, Toggle, Range1d, FreehandDrawTool
import re
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

PLOT_X_RANGE = 300
PLOT_HEIGHT = 600

def get_base_color_map():
    base_color_map = {'A': '#d6f5d6', 'C': '#ccccff', 'T': '#ffcccc', 'G': '#ffedcc', 'U': '#ffcccc', 'N': '#fafafe', 'M': '#000000', 'R': '#000000', 'W': '#000000', 'S': '#000000', 'Y': '#000000', 'K': '#000000', 'V': '#000000', 'H': '#000000', 'D': '#000000', 'B': '#000000'}
    return base_color_map

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
def create_figure(args, plot_mode):
    p_defualt = None
    if plot_mode == 0:
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
    elif plot_mode == 1:
        tools_to_show = 'hover,box_zoom,pan,save,wheel_zoom,reset,zoom_in,zoom_out'
        p_default = figure(output_backend="webgl",
                    sizing_mode="stretch_both",
                   # sizing_mode="scale_width",
                   # sizing_mode="scale_height",
                   # height=PLOT_HEIGHT,
                   x_range=(0, PLOT_X_RANGE),
                   tools=tools_to_show)
        # tooltips=tool_tips)
        # p.yaxis.visible = False
        p_default.select(dict(type=WheelZoomTool)).maintain_focus = False
        p_default.toolbar.active_scroll = p_default.select_one(WheelZoomTool)
        p_default.toolbar.logo = None
    return p_default
def scale_signal(y, sig_scale, scale_params):
    if sig_scale == "medmad":
        arr = np.ma.array(y).compressed()
        read_median = np.median(arr)
        if read_median == np.nan:
            raise Exception("Error: calculated median is NaN")
        mad = np.median(np.abs(arr - read_median))
        if mad == np.nan:
            raise Exception("Error: calculated mad is NaN")
        read_mad = mad * 1.4826
        if read_mad < 1.0:
            read_mad = 1.0
        y = (y - read_mad) / read_mad
    elif sig_scale == "znorm":
        # zsig = sklearn.preprocessing.scale(y, axis=0, with_mean=True, with_std=True, copy=True)
        # Calculate the z-score from scratch
        y = (y - np.mean(y)) / np.std(y)
    elif sig_scale == "scaledpA":
        y = (y - scale_params["sh"]) / scale_params["sc"]
    elif not sig_scale == "":
        raise Exception("Error: given --sig_scale method: {} is not supported".format(sig_scale))
    return y

profile_dic_base_shift = {
        "kmer_model_dna_r9.4.1_450bps_5_mer": [-2, -2],
        "kmer_model_dna_r9.4.1_450bps_6_mer": [-2, -3],
        "kmer_model_rna_r9.4.1_70bps_5_mer": [-3, -1],
        "kmer_model_dna_r10.4.1_e8.2_400bps_9_mer": [-6, -2],
        "guppy_dna_r9.4.1_450bps_fast_prom": [0, 0],
        "guppy_dna_r9.4.1_450bps_hac_prom": [0, 0],
        "guppy_dna_r9.4.1_450bps_sup_prom": [0, 0],
        "guppy_dna_r10.4.1_e8.2_400bps_fast": [0, 0],
        "guppy_dna_r10.4.1_e8.2_400bps_hac": [0, 0],
        "guppy_dna_r10.4.1_e8.2_400bps_sup": [0, 0]}
def list_profiles_base_shift():
    # print(profile_dic)
    print("{}\t{}\t{}".format("name", "base_shift_forward", "base_shift_reverse"))
    for profile in profile_dic_base_shift:
        print("{}\t{}\t{}".format(profile, profile_dic_base_shift[profile][0], profile_dic_base_shift[profile][1]))
    print("If the profile you wanted is not listed here, please refer calculate_offsets.md on github to learn how to generate base shift values for your new data.")
def search_for_profile_base_shift(profile):
    if profile in profile_dic_base_shift:
        return profile_dic_base_shift[profile]
    else:
        raise Exception("Error: specified profile ({}) is not found. Please run reform with -k 1 -s 0. Then run calculate_offsets.py and rerun reform with the recommended kmer_length and sig_move_offset.".format(profile))
