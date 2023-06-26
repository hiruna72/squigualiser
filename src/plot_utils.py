"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import BoxAnnotation, HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, Segment, Toggle, Range1d, FreehandDrawTool
import re
import numpy as np

PLOT_X_RANGE = 300
PLOT_HEIGHT = 600

#base_shift related
KMER_LENGTH = 9
SIG_MOVE_OFFSET = 0
BASE_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3}

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
def scale_signal(y, sig_scale):
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
    elif not sig_scale == "":
        raise Exception("Error: given --sig_scale method: {} is not supported".format(sig_scale))
    return y

def calculate_offset_values(moves, sequence, raw_signal, kmer_length, sig_move_offset):
    print("sig_move_offset: {}".format(sig_move_offset))
    print("kmer_length: {}".format(kmer_length))
    len_seq = len(sequence)
    test_array = []
    for offset in range(0, kmer_length):
        freq = [[], [], [], []]
        start_raw = 0
        for j in range(0, sig_move_offset):
            start_raw += int(moves[j])
        for i in range(0, len_seq-kmer_length + 1 - sig_move_offset):
            end_raw = start_raw + int(moves[i + sig_move_offset])
            value = raw_signal[start_raw: end_raw]
            start_raw = end_raw
            freq[BASE_INDEX[sequence[i + offset]]].append(np.median(value))
        test_array.append(freq)
    return test_array

def clean_signal(y, fasta_seq, moves):
    new_moves = []
    new_y = []
    new_fasta_seq = ""
    signal_index = 0
    base_index = 0
    print(len(fasta_seq))
    print(len(moves))
    print(fasta_seq)
    print(moves)
    for i in moves:
        if base_index >= len(fasta_seq):
            break
        if 'D' in i:
            i = re.sub('D', '', i)
            count_bases = int(i)
            base_index += count_bases
            print(base_index)

        elif 'I' in i:
            i = re.sub('I', '', i)
            eat_signal = int(i)
            signal_index += eat_signal
        else:
            eat_signal = int(i)
            count_bases = 1
            new_moves.append(i)
            for j in y[signal_index:signal_index+eat_signal]:
                new_y.append(j)
            new_fasta_seq += fasta_seq[base_index]
            base_index += count_bases
            signal_index += eat_signal
    print(102)
    return new_y, new_fasta_seq, new_moves
def calculate_base_shift(y, fasta_seq, moves):
    base_shift = 0
    kmer_length = KMER_LENGTH
    sig_move_offset = SIG_MOVE_OFFSET
    y, fasta_seq, moves = clean_signal(y, fasta_seq, moves)
    test_array = calculate_offset_values(moves, fasta_seq, y, kmer_length, sig_move_offset)
    start_offset = 0
    end_offset = kmer_length
    offset_dist = []
    for offset in range(start_offset, end_offset):
        max_mean = -1
        min_mean = 10000
        for base in test_array[offset]:
            mean = np.median(base)
            if mean < min_mean:
                min_mean = mean
            if mean > max_mean:
                max_mean = mean
        offset_dist.append(max_mean-min_mean)

    max_dist = -1
    max_offset = 0
    for offset in range(start_offset, end_offset):
        print(offset)
        print(offset_dist[offset])
        if offset_dist[offset] > max_dist:
            max_offset = offset
            max_dist = offset_dist[offset]
    print("max_offset: {}".format(max_offset))
    base_shift = -1*max_offset + sig_move_offset
    return base_shift
