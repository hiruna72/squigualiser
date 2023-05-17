"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
import numpy as np
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import BoxAnnotation, HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, Segment, Toggle, Range1d
from bokeh.layouts import row
from bokeh.colors import RGB
import pyslow5
import copy
import argparse
import re
import logging
from readpaf import parse_paf
from sys import stdin
from pyfaidx import Fasta
from pyfastx import Fastq
import os
import pysam

# ref_start is always 1based closed
# ref_end is always 1based closed
# start_kmer is always 0based closed
# end_kmer is always 0based open

FIXED_BASE_WIDTH = 10
FIXED_INSERTION_WIDTH = 10
BASE_LIMIT = 1000
SIG_PLOT_LENGTH = 20000
DEFAULT_STRIDE = 5
PLOT_X_RANGE = 750
PLOT_BOUNT_X_STRIDE = 1000
PLOT_BOUNT_Y_STRIDE = 500
PLOT_HEIGHT = 300
PLOT_BASE_SHIFT = 0

DEFAULT_NUM_BED_COLS = 3
BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK = range(10)
READ_ID, LEN_RAW_SIGNAL, START_RAW, END_RAW, STRAND, SEQUENCE_ID, LEN_KMER, START_KMER, END_KMER, MATCHES, LEN_KMER, MAPQ = range(12)
SI_START_RAW, SI_END_RAW, SI_START_KMER, SI_END_KMER = range(4)

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
        for i in moves:
            if count_bases == abs(ref_region_start_diff):
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

        moves = moves[count_moves:]
        x = signal_tuple[0][:-eat_signal]
        x_real = signal_tuple[1][eat_signal:]
        y = signal_tuple[2][eat_signal:]
        signal_tuple = (x, x_real, y)

        sig_algn_data['ss'] = moves
    return signal_tuple, region_tuple, sig_algn_data, fasta_seq
def plot_function(p, read_id, signal_tuple, sig_algn_data, fasta_sequence, base_limit, draw_data):
    x = signal_tuple[0]
    x_real = signal_tuple[1]
    y = signal_tuple[2]

    y_min = np.amin(y)
    y_max = np.amax(y)
    # label_position = np.median(y)
    # label_position = np.percentile(y, 75)  # Q3
    label_position = np.percentile(y, 98)

    base_color_map = {'A': '#d6f5d6', 'C': '#ccccff', 'T': '#ffcccc', 'G': '#ffedcc', 'U': '#ffcccc', 'N': '#fafafe'}
    base_x = []
    base_y = []
    base_label = []
    base_label_colors = []
    sample_label_colors = []
    location_plot = 0
    initial_location = location_plot

    base_shift_seq = 'N' * abs(draw_data['base_shift'])
    if draw_data["base_shift"] > 0:
        fasta_sequence = base_shift_seq + fasta_sequence[:-1*draw_data["base_shift"]]
    else:
        fasta_sequence = fasta_sequence[abs(draw_data['base_shift']):] + base_shift_seq

    # draw moves
    moves = sig_algn_data["ss"]
    base_box_details = {'left': [], 'right': [], 'fill_color': []}
    base_index = sig_algn_data["start_kmer"]
    num_Is = 0
    num_Ds = 0
    line_segment_x = []
    flag_base_index_bound = 0

    for i in moves:
        previous_location = location_plot
        if 'D' in i:
            i = re.sub('D', '', i)
            n_samples = int(i)

            prev_loc = previous_location
            for j in range(0, n_samples):
                base = fasta_sequence[base_index]

                base_box_details['left'].append(prev_loc)
                base_box_details['right'].append(prev_loc + draw_data["fixed_base_width"])
                base_box_details['fill_color'].append('white')

                base_x.append(prev_loc)
                base_y.append(label_position)
                label = str(base) + "\t" + str(base_index + 1)
                base_label.append(label)
                base_label_colors.append('red')

                prev_loc += draw_data["fixed_base_width"]
                base_index += 1
                num_Ds += 1
                if base_index - sig_algn_data["start_kmer"] == base_limit:
                    flag_base_index_bound = 1
                    break

            if flag_base_index_bound == 1:
                break
            location_plot = prev_loc

            x = x + list(range(x[-1] + 1, x[-1] + 1 + n_samples * draw_data["fixed_base_width"]))
            y_add = np.concatenate((y[:previous_location], [np.nan] * n_samples * draw_data["fixed_base_width"]), axis=0)
            y = np.concatenate((y_add, y[previous_location:]), axis=0)
            x_add = np.concatenate((x_real[:previous_location], [x_real[previous_location]] * n_samples * draw_data["fixed_base_width"]), axis=0)
            x_real = np.concatenate((x_add, x_real[previous_location:]), axis=0)
            for j in range(0, n_samples * draw_data["fixed_base_width"]):
                sample_label_colors.append('white')

        elif 'I' in i:
            i = re.sub('I', '', i)
            n_samples = int(i)
            num_Is += n_samples
            location_plot += n_samples
            line_segment_x.append(location_plot)
            for j in range(0, n_samples):
                sample_label_colors.append('purple')

            base_x.append(previous_location)
            base_y.append(label_position/2)
            label = str(n_samples)
            base_label.append(label)
            base_label_colors.append('purple')
        else:
            n_samples = int(i)
            location_plot += n_samples
            base = fasta_sequence[base_index]
            base_box_details['left'].append(previous_location)
            base_box_details['right'].append(location_plot)
            base_box_details['fill_color'].append(base_color_map[base])

            line_segment_x.append(location_plot)

            base_x.append(previous_location)
            base_y.append(label_position)
            label = str(base) + "\t" + str(base_index + 1)
            base_label.append(label)
            base_label_colors.append('black')
            for j in range(0, n_samples):
                sample_label_colors.append('red')
            base_index += 1

        if base_index - sig_algn_data["start_kmer"] == base_limit:
            break
        if location_plot - initial_location > draw_data["sig_plot_limit"]:
            break

    line_segment_source = ColumnDataSource(dict(x=line_segment_x, x1=line_segment_x, y=[y_min]*len(line_segment_x), y1=[y_max]*len(line_segment_x)))
    glyph = Segment(x0="x", y0="y", x1="x1", y1="y1", line_color="saddlebrown", line_width=1)

    base_annotation = ColumnDataSource(data=dict(base_x=base_x, base_y=base_y, base_label=base_label, colors=base_label_colors))

    base_annotation_labels = LabelSet(x='base_x', y='base_y', text='base_label', x_offset=5, y_offset=5, source=base_annotation, text_font_size="9pt", text_color='colors')

    toggle_bases = Toggle(label="base", button_type="primary", active=True, height=30, width=60)
    toggle_bases.js_link('active', base_annotation_labels, 'visible')

    source = ColumnDataSource(data=dict(x=x[:location_plot], y=y[:location_plot], x_real=x_real[:location_plot]))
    p.quad(top=y_max, bottom=y_min, left=base_box_details['left'], right=base_box_details['right'], color=base_box_details['fill_color'])
    p.add_glyph(line_segment_source, glyph)
    p.add_layout(base_annotation_labels)

    p.line('x', 'y', line_width=2, source=source)
    # add a circle renderer with a size, color, and alpha
    sample_labels = p.circle(x[:location_plot], y[:location_plot], size=draw_data["point_size"], color=sample_label_colors, alpha=0.5)
    toggle_samples = Toggle(label="sample", button_type="danger", active=True, height=30, width=60)
    toggle_samples.js_link('active', sample_labels, 'visible')

    # show the tooltip
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [("x", "@x_real"), ("y", "$y")]
    hover.mode = 'mouse'

    indt = "\t\t\t\t\t\t\t\t"

    if sig_algn_data["plot_sig_ref_flag"] == 0:
        if sig_algn_data["data_is_rna"] == 1:
            plot_title = f'base_shift: {draw_data["base_shift"]}{indt}{sig_algn_data["tag_name"]}[{base_index}-{sig_algn_data["ref_start"]}]{indt}signal: [{int(x_real[0])}-{int(x_real[location_plot - 1])}]{indt}deletions(bases): {num_Ds} insertions(samples): {num_Is}{indt}{read_id}'
        else:
            plot_title = f'base_shift: {draw_data["base_shift"]}{indt}{sig_algn_data["tag_name"]}[{sig_algn_data["ref_start"]}-{base_index}]{indt}signal: [{int(x_real[0])}-{int(x_real[location_plot - 1])}]{indt}deletions(bases): {num_Ds} insertions(samples): {num_Is}{indt}{read_id}'
    else:
        if sig_algn_data["data_is_rna"] == 1:
            plot_title = f'base_shift: {draw_data["base_shift"]}{indt}{sig_algn_data["tag_name"]}[{sig_algn_data["ref_end"]}-{sig_algn_data["ref_end"] - base_index + 1}]{indt}signal: [{int(x_real[0])}-{int(x_real[location_plot - 1])}]{indt}deletions(bases): {num_Ds} insertions(samples): {num_Is}{indt}{read_id}'
        else:
            plot_title = f'base_shift: {draw_data["base_shift"]}{indt}{sig_algn_data["tag_name"]}[{sig_algn_data["ref_start"]}-{sig_algn_data["ref_start"] + base_index - 1}]{indt}signal: [{int(x_real[0])}-{int(x_real[location_plot - 1])}]{indt}deletions(bases): {num_Ds} insertions(samples): {num_Is}{indt}{read_id}'
    p.title = plot_title

    layout_ = p, row(toggle_bases, toggle_samples)
    return layout_
def plot_function_fixed_width(p, read_id, signal_tuple, sig_algn_data, fasta_sequence, base_limit, draw_data):
    x = signal_tuple[0]
    x_real = signal_tuple[1]
    y = signal_tuple[2]

    # label_position = np.median(y)
    # label_position = np.percentile(y, 75)  # Q3
    label_position = np.percentile(y, 98)
    y_min = np.amin(y)
    y_max = np.amax(y)

    base_color_map = {'A': '#d6f5d6', 'C': '#ccccff', 'T': '#ffcccc', 'G': '#ffedcc', 'U': '#ffcccc', 'N': '#fafafe'}
    base_x = []
    base_y = []
    base_label = []
    base_label_colors = []
    sample_label_colors = []
    location_plot = 0
    initial_location = location_plot

    x_coordinate = 0
    initial_x_coordinate = x_coordinate

    base_shift_seq = 'N' * abs(draw_data['base_shift'])
    if draw_data["base_shift"] > 0:
        fasta_sequence = base_shift_seq + fasta_sequence[:-1*draw_data["base_shift"]]
    else:
        fasta_sequence = fasta_sequence[abs(draw_data['base_shift']):] + base_shift_seq

    # draw moves
    moves = sig_algn_data["ss"]
    base_index = sig_algn_data["start_kmer"]
    num_Is = 0
    num_Ds = 0
    fixed_width_x = [0.0]
    line_segment_x = []
    base_box_details = {'left': [], 'right': [], 'fill_color': []}

    num_samples_in_insertion = 0
    flag_base_index_bound = 0
    for i in moves:
        previous_location = location_plot
        previous_x_coordinate = x_coordinate
        if 'D' in i:
            i = re.sub('D', '', i)
            n_samples = int(i)

            prev_loc = previous_location
            prev_x_cord = previous_x_coordinate
            for j in range(0, n_samples):
                base = fasta_sequence[base_index]
                base_box_details['left'].append(prev_loc)
                base_box_details['right'].append(prev_loc + draw_data["fixed_base_width"])
                base_box_details['fill_color'].append('white')

                base_x.append(prev_loc)
                base_y.append(label_position)
                label = str(base) + "\t" + str(base_index + 1)
                base_label.append(label)
                base_label_colors.append('red')

                prev_loc += draw_data["fixed_base_width"]
                prev_x_cord += draw_data["fixed_base_width"]
                base_index += 1
                num_Ds += 1
                if base_index - sig_algn_data["start_kmer"] == base_limit:
                    flag_base_index_bound = 1
                    break

            if flag_base_index_bound == 1:
                break
            location_plot = prev_loc
            x_coordinate = prev_x_cord

            fixed_width_x = fixed_width_x + list(range(int(fixed_width_x[-1]) + 1, int(fixed_width_x[-1]) + 1 + n_samples * draw_data["fixed_base_width"]))
            y_add = np.concatenate((y[:previous_x_coordinate], [np.nan] * n_samples * draw_data["fixed_base_width"]), axis=0)
            y = np.concatenate((y_add, y[previous_x_coordinate:]), axis=0)
            x_add = np.concatenate((x_real[:previous_x_coordinate], [x_real[previous_x_coordinate]] * n_samples * draw_data["fixed_base_width"]), axis=0)
            x_real = np.concatenate((x_add, x_real[previous_x_coordinate:]), axis=0)

            for j in range(0, n_samples * draw_data["fixed_base_width"]):
                sample_label_colors.append('white')

        elif 'I' in i:
            i = re.sub('I', '', i)
            n_samples = int(i)
            # prev_x_value = fixed_width_x[-1]
            # for s in range(0, n_samples-1):
            #     fixed_width_x.append(fixed_width_x[-1] + draw_data["fixed_base_width"]/n_samples)
            # fixed_width_x.append(prev_x_value + draw_data["fixed_base_width"])
            num_samples_in_insertion = n_samples
            num_Is += n_samples
            # location_plot += draw_data["fixed_base_width"]
            # x_coordinate += n_samples
            # line_segment_x.append(location_plot)
            for j in range(0, n_samples):
                sample_label_colors.append('purple')

        else:
            n_samples = int(i) + num_samples_in_insertion
            prev_x_value = fixed_width_x[-1]
            for s in range(0, n_samples - 1):
                fixed_width_x.append(fixed_width_x[-1] + draw_data["fixed_base_width"] / n_samples)
            fixed_width_x.append(prev_x_value + draw_data["fixed_base_width"])

            location_plot += draw_data["fixed_base_width"]
            x_coordinate += n_samples

            base = fasta_sequence[base_index]
            base_box_details['left'].append(previous_location)
            base_box_details['right'].append(location_plot)
            base_box_details['fill_color'].append(base_color_map[base])
            line_segment_x.append(location_plot)

            if num_samples_in_insertion > 0:
                base_x.append(previous_location)
                base_y.append(label_position/2)
                label = str(num_samples_in_insertion)
                base_label.append(label)
                base_label_colors.append('purple')

            base_x.append(previous_location)
            base_y.append(label_position)
            label = str(base) + "\t" + str(base_index + 1)
            base_label.append(label)
            base_label_colors.append('black')
            for j in range(0, n_samples - num_samples_in_insertion):
                sample_label_colors.append('red')
            base_index += 1
            num_samples_in_insertion = 0

        if base_index - sig_algn_data["start_kmer"] == base_limit:
            break
        if x_coordinate - initial_x_coordinate > draw_data["sig_plot_limit"]:
            break

    line_segment_source = ColumnDataSource(dict(x=line_segment_x, x1=line_segment_x, y=[y_min] * len(line_segment_x), y1=[y_max] * len(line_segment_x)))
    glyph = Segment(x0="x", y0="y", x1="x1", y1="y1", line_color="saddlebrown", line_width=1)

    base_annotation = ColumnDataSource(data=dict(base_x=base_x, base_y=base_y, base_label=base_label, colors=base_label_colors))

    base_annotation_labels = LabelSet(x='base_x', y='base_y', text='base_label', x_offset=5, y_offset=5, source=base_annotation, text_font_size="9pt", text_color='colors')

    toggle_bases = Toggle(label="base", button_type="primary", active=True, height=30, width=60)
    toggle_bases.js_link('active', base_annotation_labels, 'visible')

    fixed_width_x = fixed_width_x[1:]

    source = ColumnDataSource(data=dict(x=fixed_width_x[:x_coordinate], y=y[:x_coordinate], x_real=x_real[:x_coordinate]))
    p.quad(top=y_max, bottom=y_min, left=base_box_details['left'], right=base_box_details['right'], color=base_box_details['fill_color'])
    p.add_glyph(line_segment_source, glyph)
    p.add_layout(base_annotation_labels)

    p.line('x', 'y', line_width=2, source=source)
    # add a circle renderer with a size, color, and alpha
    sample_labels = p.circle(fixed_width_x[:x_coordinate], y[:x_coordinate], size=draw_data["point_size"], color=sample_label_colors, alpha=0.5)
    toggle_samples = Toggle(label="sample", button_type="danger", active=True, height=30, width=60)
    toggle_samples.js_link('active', sample_labels, 'visible')

    # show the tooltip
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [("x", "@x_real"), ("y", "$y")]
    hover.mode = 'mouse'

    indt = "\t\t\t\t\t\t\t\t"

    if sig_algn_data["plot_sig_ref_flag"] == 0:
        if sig_algn_data["data_is_rna"] == 1:
            plot_title = f'base_shift: {draw_data["base_shift"]}{indt}{sig_algn_data["tag_name"]}[{base_index}-{sig_algn_data["ref_start"]}]{indt}signal: [{int(x_real[0])}-{int(x_real[x_coordinate - 1])}]{indt}deletions(bases): {num_Ds} insertions(samples): {num_Is}{indt}{read_id}'
        else:
            plot_title = f'base_shift: {draw_data["base_shift"]}{indt}{sig_algn_data["tag_name"]}[{sig_algn_data["ref_start"]}-{base_index}]{indt}signal: [{int(x_real[0])}-{int(x_real[x_coordinate - 1])}]{indt}deletions(bases): {num_Ds} insertions(samples): {num_Is}{indt}{read_id}'
    else:
        if sig_algn_data["data_is_rna"] == 1:
            plot_title = f'base_shift: {draw_data["base_shift"]}{indt}{sig_algn_data["tag_name"]}[{sig_algn_data["ref_end"]}-{sig_algn_data["ref_end"] - base_index + 1}]{indt}signal: [{int(x_real[0])}-{int(x_real[x_coordinate - 1])}]{indt}deletions(bases): {num_Ds} insertions(samples): {num_Is}{indt}{read_id}'
        else:
            plot_title = f'base_shift: {draw_data["base_shift"]}{indt}{sig_algn_data["tag_name"]}[{sig_algn_data["ref_start"]}-{sig_algn_data["ref_start"] + base_index - 1}]{indt}signal: [{int(x_real[0])}-{int(x_real[x_coordinate - 1])}]{indt}deletions(bases): {num_Ds} insertions(samples): {num_Is}{indt}{read_id}'
    p.title = plot_title

    if location_plot < PLOT_X_RANGE:
        draw_data["plot_dims"]['start_x'] = 0
        draw_data["plot_dims"]['end_x'] = PLOT_X_RANGE
    else:
        draw_data["plot_dims"]['start_x'] = location_plot/2 - PLOT_X_RANGE/2
        draw_data["plot_dims"]['end_x'] = location_plot/2 + PLOT_X_RANGE/2
    draw_data["plot_dims"]['start_y'] = y_min - (y_max-y_min)/2
    draw_data["plot_dims"]['end_y'] = y_max + (y_max-y_min)/2
    draw_data["plot_dims"]['bound_start_x'] = 0 - PLOT_BOUNT_X_STRIDE
    draw_data["plot_dims"]['bound_end_x'] = location_plot + PLOT_BOUNT_X_STRIDE
    draw_data["plot_dims"]['bound_start_y'] = y_min - PLOT_BOUNT_Y_STRIDE
    draw_data["plot_dims"]['bound_end_y'] = y_max + PLOT_BOUNT_Y_STRIDE
    p.x_range = Range1d(draw_data["plot_dims"]['start_x'], draw_data["plot_dims"]['end_x'], bounds=(draw_data["plot_dims"]['bound_start_x'], draw_data["plot_dims"]['bound_end_x']))
    p.y_range = Range1d(draw_data["plot_dims"]['start_y'], draw_data["plot_dims"]['end_y'], bounds=(draw_data["plot_dims"]['bound_start_y'], draw_data["plot_dims"]['bound_end_y']))

    layout_ = p, row(toggle_bases, toggle_samples)
    return layout_

def plot_bed_annotation(p, annotation, signal_tuple):
    p.quad(top=400, bottom=300, left=0, right=10, color='red')
    return p

def run(args):
    if args.read_id != "":
        args.plot_limit = 1

    use_fasta = 0
    if args.file:
        print(f'sequence file: {args.file}')
        if args.file[-6:] == ".fastq" or args.file[-6:] == ".fq.gz" or args.file[-3:] == ".fq":
            use_fasta = 0
        elif args.file[-6:] == ".fasta" or args.file[-3:] == ".fa":
            use_fasta = 1
        else:
            raise Exception("Error: please provide the sequence file with correct extension")

    plot_sig_ref_flag = 0
    if args.sig_ref:
        plot_sig_ref_flag = 1
    use_paf = 0
    if args.alignment:
        print(f'alignment file: {args.alignment}')
        alignment_extension = args.alignment[-4:]
        if alignment_extension == ".bam" or alignment_extension == ".sam":
            use_paf = 0
            plot_sig_ref_flag = 1
        elif alignment_extension == ".paf":
            use_paf = 1
            if args.sig_ref:
                raise Exception("Error: For signal to reference alignment the paf file must be bgzip compressed and tabix indexed")
        elif args.alignment[-7:] == ".paf.gz":
            use_paf = 1
            plot_sig_ref_flag = 1
            index_file = args.alignment + ".tbi"
            if not os.path.exists(index_file):
                raise Exception("Error: please provide a bgzip compressed and tabix indexed paf file to extract the regions")
        else:
            raise Exception("Error: please provide the alignment file with correct extension")

    if use_paf == 0 and use_fasta == 0:
        print("please provide a .fasta or .fa file when using SAM/BAM")


    if args.base_limit:
        base_limit = args.base_limit
    else:
        base_limit = BASE_LIMIT
    print(f'signal file: {args.slow5}')

    bed_content = []
    if args.bed:
        print(f'bed file: {args.bed}')
        bed_num_cols = DEFAULT_NUM_BED_COLS
        with open(args.bed)as f:
            for line in f:
                bed_list = line.strip().split()
                if len(bed_list) < bed_num_cols:
                    raise Exception("Error: minimum {} columns required in a bed file", bed_num_cols)
                else:
                    bed_num_cols = len(bed_list)
                bed_content.append(bed_list)
        print(bed_content)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # open signal file
    s5 = pyslow5.Open(args.slow5, 'r')
    num_plots = 0
    indt = "\t\t\t\t\t\t\t\t"
    draw_data = {}
    draw_data["point_size"] = args.point_size
    draw_data["sig_plot_limit"] = args.sig_plot_limit
    draw_data["fixed_base_width"] = args.base_width
    draw_data["base_shift"] = args.base_shift
    draw_data["plot_dims"] = {}
    draw_data["bed_content"] = bed_content

    y_axis_label = "signal value (raw)"
    if args.no_pa:
        y_axis_label = "signal value (pA)"
    tools_to_show = 'hover,box_zoom,pan,reset,save,wheel_zoom,zoom_in,zoom_out'
    p = figure(x_axis_label='signal index',
               y_axis_label=y_axis_label,
               sizing_mode="stretch_width",
               height=PLOT_HEIGHT,
               output_backend="webgl",
               x_range=(0, PLOT_X_RANGE),
               tools=tools_to_show,
               toolbar_location="below")
    # tooltips=tool_tips)
    p.toolbar.active_scroll = p.select_one(WheelZoomTool)
    p.toolbar.logo = None

    if use_paf == 1 and plot_sig_ref_flag == 0:
        print("Info: Signal to read method using PAF ...")
        with open(args.alignment, "r") as handle:
            if use_fasta:
                sequence_reads = Fasta(args.file)
            else:
                sequence_reads = Fastq(args.file)
            for paf_record in parse_paf(handle):
                if paf_record.query_name != paf_record.target_name:
                    raise Exception("Error: this paf file is a signal to reference mapping. Please provide the argument --sig_ref ")
                read_id = paf_record.query_name
                if args.read_id != "" and read_id != args.read_id:
                    continue
                if read_id not in set(sequence_reads.keys()):
                    raise Exception("Error: read_id {} is not found in {}".format(read_id, args.file))

                data_is_rna = 0
                if paf_record.target_start > paf_record.target_end:  # if RNA start_kmer>end_kmer in paf
                    data_is_rna = 1
                    if not args.rna:
                        print("Info: data is detected as RNA")
                        raise Exception("Error: data is not specified as RNA. Please provide the argument --rna ")

                fasta_seq = ""
                if use_fasta:
                    fasta_seq = sequence_reads[read_id][:].seq
                else:
                    fasta_seq = sequence_reads[read_id].seq
                    if len(fasta_seq) < paf_record.target_length:
                        print(
                            "Error: Sequence lengths mismatch. If {} is a multi-line fastq file convert it to a 4-line fastq using seqtk.".format(
                                args.file))
                        exit(1)

                ref_start = -1
                ref_end = -1
                if data_is_rna == 1:
                    fasta_seq = fasta_seq[paf_record.target_end:]
                    ref_start = paf_record.target_end + 1
                else:
                    fasta_seq = fasta_seq[paf_record.target_start:]
                    ref_start = paf_record.target_start + 1

                seq_len = len(fasta_seq)
                if seq_len < base_limit:
                    base_limit = seq_len
                ref_end = base_limit

                if args.region != "":
                    pattern = re.compile("^[0-9]+\-[0-9]+")
                    if not pattern.match(args.region):
                        raise Exception("Error: region provided is not in correct format")
                    ref_start = int(args.region.split("-")[0])
                    ref_end = int(args.region.split("-")[1])

                    if ref_start < 1:
                        raise Exception("Error: region start coordinate ({}) must be positive".format(ref_start))

                    if ref_end < 1:
                        raise Exception("Error: region end coordinate ({}) must be positive".format(ref_end))

                    if data_is_rna == 1 and paf_record.target_start < ref_end:
                        ref_end = paf_record.target_start
                    elif data_is_rna == 0 and paf_record.target_end < ref_end:
                        ref_end = paf_record.target_end

                    # print("ref_start: " + str(ref_start))
                    # print("ref_end: " + str(ref_end))

                    if (ref_end - ref_start + 1) < base_limit:
                        base_limit = ref_end - ref_start + 1

                print("plot region: {}-{}\tread_id: {}".format(ref_start, ref_end, read_id))

                output_file_name = args.output_dir + "/" + read_id + "_" + args.tag_name + ".html"

                x = []
                x_real = []
                y = []
                read = s5.get_read(read_id, pA=args.no_pa, aux=["read_number", "start_mux"])
                if read is not None:
                    start_index = paf_record.query_start
                    end_index = read['len_raw_signal']
                    x = list(range(1, end_index - start_index + 1))
                    x_real = list(range(start_index + 1, end_index + 1))  # 1based
                    y = read['signal'][start_index:end_index]
                if args.sig_scale == "medmad":
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
                    args.tag_name += " scale:medmad"
                elif args.sig_scale == "znorm":
                    # zsig = sklearn.preprocessing.scale(y, axis=0, with_mean=True, with_std=True, copy=True)
                    # Calculate the z-score from scratch
                    y = (y - np.mean(y)) / np.std(y)
                    args.tag_name += " scale:znorm"
                elif not args.sig_scale == "":
                    raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))

                strand_dir = "(DNA 5'->3')"

                # if data_is_rna == 1 and args.reverse_signal:
                #     x_real.reverse()
                #     y = np.flip(y)
                #     moves.reverse()
                #     strand_dir = "(RNA 5'->3')"

                if data_is_rna == 1:
                    fasta_seq = fasta_seq[::-1]
                    strand_dir = "(RNA 5'->3')"

                signal_tuple = (x, x_real, y)
                region_tuple = (ref_start, ref_end, 0, seq_len)

                sig_algn_dic = {}
                sig_algn_dic['start_kmer'] = ref_start - 1
                sig_algn_dic['ref_start'] = ref_start
                sig_algn_dic['ref_end'] = ref_end
                sig_algn_dic['use_paf'] = use_paf
                sig_algn_dic['plot_sig_ref_flag'] = plot_sig_ref_flag
                sig_algn_dic['data_is_rna'] = data_is_rna
                if args.fixed_width:
                    sig_algn_dic['tag_name'] = args.tag_name + indt + "fixed_width: " + str(args.base_width) + indt + strand_dir + indt + "region: "
                else:
                    sig_algn_dic['tag_name'] = args.tag_name + indt + strand_dir + indt + "region: "

                moves_string = paf_record.tags['ss'][2]
                moves_string = re.sub('D', 'D,', moves_string)
                moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
                moves = re.split(r',+', moves_string)
                sig_algn_dic['ss'] = moves

                signal_tuple, region_tuple, sig_algn_dic, fasta_seq = adjust_before_plotting(seq_len, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)
                if args.fixed_width:
                    layout_ = plot_function_fixed_width(p, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)
                else:
                    layout_ = plot_function(p=p, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)

                output_file(output_file_name, title=read_id)
                save(layout_)
                print(f'output file: {os.path.abspath(output_file_name)}')

                num_plots += 1
                if num_plots == args.plot_limit:
                    break
    elif use_paf == 0 and plot_sig_ref_flag == 1: # using sam/bam
        print("Info: Signal to reference method using SAM/BAM ...")
        fasta_reads = Fasta(args.file)
        if args.region != "":
            # check if there exists a .bam.bai
            index_file = args.alignment + ".bai"
            if not os.path.exists(index_file):
                raise Exception("Error: please provide a bam file that is sorted and indexed to extract the regions")

            args_region = re.sub(',', '', args.region)
            # print(args_region)
            # pattern = re.compile("^[a-z]+[0-9]+\:[0-9]+\-[0-9]+")
            pattern = re.compile("^.*\:[0-9]+\-[0-9]+")
            if not pattern.match(args_region):
                raise Exception("Error: region provided is not in correct format")
            args_ref_name = args_region.split(":")[0]
            args_ref_start = int(args_region.split(":")[1].split("-")[0])
            args_ref_end = int(args_region.split(":")[1].split("-")[1])

            samfile = pysam.AlignmentFile(args.alignment, mode='rb')
        else:
            args_ref_name = None
            args_ref_start = None
            args_ref_end = None
            samfile = pysam.AlignmentFile(args.alignment, mode='r')

        for sam_record in samfile.fetch(contig=args_ref_name, start=args_ref_start, stop=args_ref_end):
            if args_ref_name is not None and args_ref_name != sam_record.reference_name:
                raise Exception("Error: sam record's reference name [" + sam_record.reference_name + "] and the name specified are different [" + args_ref_name + "]")
            read_id = sam_record.query_name
            if sam_record.is_supplementary or sam_record.is_unmapped or sam_record.is_secondary:
                continue
            if sam_record.is_reverse and args.no_reverse:
                continue
            if not sam_record.is_reverse and args.reverse_only:
                continue
            if args.read_id != "" and read_id != args.read_id:
                continue

            ref_seq_len = 0
            data_is_rna = 0
            start_index = -1
            if sam_record.has_tag("si"):
                si_tag = sam_record.get_tag("si").split(',')
                start_index = int(si_tag[SI_START_RAW])
                end_index = int(si_tag[SI_END_RAW])
                ref_seq_len = int(si_tag[SI_END_KMER]) - int(si_tag[SI_START_KMER])

                if int(si_tag[SI_START_KMER]) > int(si_tag[SI_END_KMER]):  # if RNA start_kmer>end_kmer in paf
                    data_is_rna = 1
                    if not args.rna:
                        print("Info: data is detected as RNA")
                        raise Exception("Error: data is not specified as RNA. Please provide the argument --rna ")
                    ref_seq_len = int(si_tag[SI_START_KMER]) - int(si_tag[SI_END_KMER])

            else:
                raise Exception("Error: sam record does not have a 'si' tag.")
            # print("ref_seq_len: " + str(ref_seq_len))

            if ref_seq_len < BASE_LIMIT:
                base_limit = ref_seq_len
            else:
                base_limit = BASE_LIMIT
            sam_record_reference_end = sam_record.reference_start + ref_seq_len #1based closed

            ref_name = sam_record.reference_name
            ref_start = sam_record.reference_start + 1
            ref_end = ref_start + base_limit - 1 #ref_end is 1based closed
            if args.region != "":
                if args_ref_start > ref_start:
                    ref_start = args_ref_start
                    if (ref_start + base_limit - 1) < sam_record_reference_end:
                        ref_end = ref_start + base_limit - 1
                    else:
                        ref_end = sam_record_reference_end
                if args_ref_end < ref_end:
                    ref_end = args_ref_end
            if ref_end < ref_start:
                print("Warning: a corner case has hit because  the kmer_length used is larger than 1. This alignment will be skipped")
                continue
            base_limit = ref_end - ref_start + 1
            # print("ref_start: {}".format(ref_start))
            # print("ref_end: {}".format(ref_end))
            # print("ref_seq_len: {}".format(ref_seq_len))
            # print("{}".format(sam_record.cigarstring))
            # print("base_limit: {}".format(base_limit))

            if data_is_rna == 1:
                print("plot (RNA 5'->3') region: {}:{}-{}\tread_id: {}".format(ref_name, ref_end, ref_start, read_id))
                fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
            else:
                if sam_record.is_reverse:
                    print("plot (DNA 5'->3' -) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                else:
                    print("plot (DNA 5'->3' +) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
            output_file_name = args.output_dir + "/" + read_id + "_" + args.tag_name + ".html"

            x = []
            x_real = []
            y = []

            read = s5.get_read(read_id, pA=args.no_pa, aux=["read_number", "start_mux"])
            if read is not None:
                # print("read_id:", read['read_id'])
                # print("len_raw_signal:", read['len_raw_signal'])
                # end_index = read['len_raw_signal']
                x = list(range(1, end_index - start_index + 1))
                x_real = list(range(start_index+1, end_index+1))             # 1based
                y = read['signal'][start_index:end_index]
            if args.sig_scale == "medmad":
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
                args.tag_name += " scale:medmad"
            elif args.sig_scale == "znorm":
                # zsig = sklearn.preprocessing.scale(y, axis=0, with_mean=True, with_std=True, copy=True)
                # Calculate the z-score from scratch
                y = (y - np.mean(y)) / np.std(y)
                args.tag_name += " scale:znorm"
            elif not args.sig_scale == "":
                raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))

            moves_string = sam_record.get_tag("ss")
            moves_string = re.sub('D', 'D,', moves_string)
            moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
            moves = re.split(r',+', moves_string)

            # if data_is_rna == 0 and args.reverse_signal:
            #     raise Exception("Error: the signal will be reversed only in RNA plots")

            if sam_record.is_reverse and data_is_rna == 1:
                raise Exception("Error: the signal is  always sequenced from 3` to 5`. Hence, cannot have reversed mapped reads?")

            # if data_is_rna == 1 and args.reverse_signal:
            #     x_real.reverse()
            #     y = np.flip(y)
            #     moves.reverse()
            #     strand_dir = "(RNA 5'->3')"
            if data_is_rna == 0:
                strand_dir = "(DNA +)"
                if sam_record.is_reverse:
                    strand_dir = "(DNA -)"
                    x_real.reverse()
                    y = np.flip(y)
                    moves.reverse()

            if data_is_rna == 1:
                strand_dir = "(RNA 5'->3')"
                fasta_seq = fasta_seq[::-1]
                if sam_record.is_reverse:
                    raise Exception("Error: data is rna and sam record is reverse mapped. This is not implemented yet. Please report")

            signal_tuple = (x, x_real, y)
            region_tuple = (ref_start, ref_end, sam_record.reference_start, sam_record.reference_start+ref_seq_len)

            sig_algn_dic = {}
            sig_algn_dic['start_kmer'] = 0
            sig_algn_dic['ref_start'] = ref_start
            sig_algn_dic['ref_end'] = ref_end
            sig_algn_dic['plot_sig_ref_flag'] = plot_sig_ref_flag
            sig_algn_dic['data_is_rna'] = data_is_rna
            if args.fixed_width:
                sig_algn_dic['tag_name'] = args.tag_name + indt + "fixed_width: " + str(args.base_width) + indt + strand_dir + indt + "region: " + ref_name + ":"
            else:
                sig_algn_dic['tag_name'] = args.tag_name + indt + strand_dir + indt + "region: " + ref_name + ":"
            sig_algn_dic['ss'] = moves
            # print(len(moves))
            # print(fasta_seq)
            signal_tuple, region_tuple, sig_algn_dic, fasta_seq = adjust_before_plotting(ref_seq_len, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)
            # print(len(sig_algn_dic['ss']))

            if args.fixed_width:
                layout_ = plot_function_fixed_width(p, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)
            else:
                layout_ = plot_function(p=p, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)

            output_file(output_file_name, title=read_id)
            save(layout_)
            print(f'output file: {os.path.abspath(output_file_name)}')
                        
            num_plots += 1
            if num_plots == args.plot_limit:
                break
    elif use_paf == 1 and plot_sig_ref_flag == 1:
        print("Info: Signal to reference method using PAF ...")
        fasta_reads = Fasta(args.file)
        tbxfile = pysam.TabixFile(args.alignment)
        if args.region != "":
            args_region = re.sub(',', '', args.region)
            # print(args_region)
            # pattern = re.compile("^[a-z]+[0-9]+\:[0-9]+\-[0-9]+")
            pattern = re.compile("^.*\:[0-9]+\-[0-9]+")
            if not pattern.match(args_region):
                raise Exception("Error: region provided is not in correct format")
            args_ref_name = args_region.split(":")[0]
            args_ref_start = int(args_region.split(":")[1].split("-")[0])
            args_ref_end = int(args_region.split(":")[1].split("-")[1])
        else:
            args_ref_name = None
            args_ref_start = None
            args_ref_end = None

        for paf_record in tbxfile.fetch(args_ref_name, args_ref_start, args_ref_end, parser=pysam.asTuple()):
            if paf_record[READ_ID] == paf_record[SEQUENCE_ID]:
                raise Exception("Error: this paf file is a signal to read mapping.")
            if args_ref_name is not None and args_ref_name != paf_record[SEQUENCE_ID]:
                raise Exception("Error: sam record's reference name [" + paf_record[SEQUENCE_ID] + "] and the name specified are different [" + args_ref_name + "]")
            read_id = paf_record[READ_ID]
            # if read_id != "285802f0-8f4d-4f03-8d11-ef8a395576e4":
            #     continue
            if args.read_id != "" and read_id != args.read_id:
                continue
            if paf_record[STRAND] == "-" and args.no_reverse:
                continue
            if paf_record[STRAND] == "+" and args.reverse_only:
                continue

            data_is_rna = 0
            start_index = int(paf_record[START_RAW])
            end_index = int(paf_record[END_RAW])
            ref_seq_len = int(paf_record[END_KMER]) - int(paf_record[START_KMER])
            reference_start = int(paf_record[START_KMER])
            if int(paf_record[START_KMER]) > int(paf_record[END_KMER]):  # if RNA start_kmer>end_kmer in paf
                data_is_rna = 1
                if not args.rna:
                    print("Info: data is detected as RNA")
                    raise Exception("Error: data is not specified as RNA. Please provide the argument --rna ")
                ref_seq_len = int(paf_record[START_KMER]) - int(paf_record[END_KMER])
                reference_start = int(paf_record[END_KMER])
            # print("ref_seq_len: " + str(ref_seq_len))
            if ref_seq_len < BASE_LIMIT:
                base_limit = ref_seq_len
            else:
                base_limit = BASE_LIMIT
            paf_record_reference_end = reference_start + ref_seq_len #1based closed

            ref_name = paf_record[SEQUENCE_ID]
            ref_start = reference_start + 1
            ref_end = ref_start + base_limit - 1 #ref_end is 1based closed
            if args.region != "":
                if args_ref_start > ref_start:
                    ref_start = args_ref_start
                    if (ref_start + base_limit - 1) < paf_record_reference_end:
                        ref_end = ref_start + base_limit - 1
                    else:
                        ref_end = paf_record_reference_end
                if args_ref_end < ref_end:
                    ref_end = args_ref_end
            if ref_end < ref_start:
                print("Warning: a corner case has hit because the kmer_length used is larger than 1. This alignment will be skipped")
                continue
            base_limit = ref_end - ref_start + 1
            # print("ref_start: {}".format(ref_start))
            # print("ref_end: {}".format(ref_end))
            # print("ref_seq_len: {}".format(ref_seq_len))
            # print("base_limit: {}".format(base_limit))

            record_is_reverse = 0
            if paf_record[STRAND] == "-":
                record_is_reverse = 1
            if record_is_reverse and data_is_rna == 1:
                raise Exception("Error: the signal is  always sequenced from 3` to 5`. Hence, cannot have reversed mapped reads?")
            if data_is_rna == 1:
                print("plot (RNA 5'->3') region: {}:{}-{}\tread_id: {}".format(ref_name, ref_end, ref_start, read_id))
                fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
            else:
                if record_is_reverse:
                    print("plot (DNA 5'->3' -) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                else:
                    print("plot (DNA 5'->3' +) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
            output_file_name = args.output_dir + "/" + read_id + "_" + args.tag_name + ".html"

            x = []
            x_real = []
            y = []

            read = s5.get_read(read_id, pA=args.no_pa, aux=["read_number", "start_mux"])
            if read is not None:
                # print("read_id:", read['read_id'])
                # print("len_raw_signal:", read['len_raw_signal'])
                # end_index = read['len_raw_signal']
                x = list(range(1, end_index - start_index + 1))
                x_real = list(range(start_index + 1, end_index + 1))  # 1based
                y = read['signal'][start_index:end_index]
            if args.sig_scale == "medmad":
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
                args.tag_name += " scale:medmad"
            elif args.sig_scale == "znorm":
                # zsig = sklearn.preprocessing.scale(y, axis=0, with_mean=True, with_std=True, copy=True)
                # Calculate the z-score from scratch
                y = (y - np.mean(y)) / np.std(y)
                args.tag_name += " scale:znorm"
            elif not args.sig_scale == "":
                raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))

            for i in range(12, len(paf_record)):
                tag = paf_record[i][:2]
                if tag == "ss":
                    moves_string = paf_record[i][5:]
            if moves_string == "":
                raise Exception("Error: ss tag was not found")
            moves_string = re.sub('D', 'D,', moves_string)
            moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
            moves = re.split(r',+', moves_string)

            if data_is_rna == 0:
                strand_dir = "(DNA +)"
                if record_is_reverse:
                    strand_dir = "(DNA -)"
                    x_real.reverse()
                    y = np.flip(y)
                    moves.reverse()

            if data_is_rna == 1:
                strand_dir = "(RNA 5'->3')"
                fasta_seq = fasta_seq[::-1]

            signal_tuple = (x, x_real, y)
            region_tuple = (ref_start, ref_end, reference_start, reference_start + ref_seq_len)

            sig_algn_dic = {}
            sig_algn_dic['start_kmer'] = 0
            sig_algn_dic['ref_start'] = ref_start
            sig_algn_dic['ref_end'] = ref_end
            sig_algn_dic['pa'] = args.no_pa
            sig_algn_dic['plot_sig_ref_flag'] = plot_sig_ref_flag
            sig_algn_dic['data_is_rna'] = data_is_rna
            if args.fixed_width:
                sig_algn_dic['tag_name'] = args.tag_name + indt + "fixed_width: " + str(
                    args.base_width) + indt + strand_dir + indt + "region: " + ref_name + ":"
            else:
                sig_algn_dic['tag_name'] = args.tag_name + indt + strand_dir + indt + "region: " + ref_name + ":"
            sig_algn_dic['ss'] = moves
            # print(len(moves))
            # print(fasta_seq)
            signal_tuple, region_tuple, sig_algn_dic, fasta_seq = adjust_before_plotting(ref_seq_len, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)
            # print(len(sig_algn_dic['ss']))

            if args.fixed_width:
                layout_ = plot_function_fixed_width(p, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)
            else:
                layout_ = plot_function(p=p, read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data)
                
            output_file(output_file_name, title=read_id)
            save(layout_)
            print(f'output file: {os.path.abspath(output_file_name)}')

            num_plots += 1
            if num_plots == args.plot_limit:
                break
    else:
        raise Exception("Error: You should not have ended up here. Please check your arguments")



    print("Number of plots: {}".format(num_plots))

    s5.close()

def argparser():
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    parser.add_argument('-f', '--file', required=True, help="fasta/fa/fastq/fq/fq.gz sequence file")
    parser.add_argument('-r', '--read_id', required=False, type=str, default="", help="plot the read with read_id")
    parser.add_argument('--base_limit', required=False, type=int, help="maximum number of bases to plot")
    parser.add_argument('-s', '--slow5', required=True, help="slow5 file")
    parser.add_argument('-a', '--alignment', required=True,
                        help="for read-signal alignment use PAF\nfor reference-signal alignment use SAM/BAM")
    parser.add_argument('--region', required=False, type=str, default="",
                        help="[start-end] 1-based closed interval region to plot. For SAM/BAM eg: chr1:6811428-6811467 or chr1:6,811,428-6,811,467. For PAF eg:100-200.")
    parser.add_argument('--tag_name', required=False, type=str, default="",
                        help="a tag name to easily identify the plot")
    parser.add_argument('--no_reverse', required=False, action='store_true', help="skip plotting reverse mapped reads")
    parser.add_argument('--reverse_only', required=False, action='store_true', help="only plot reverse mapped reads")
    parser.add_argument('--rna', required=False, action='store_true', help="specify for RNA reads")
    parser.add_argument('--sig_ref', required=False, action='store_true', help="plot signal to reference mapping")
    parser.add_argument('--fixed_width', required=False, action='store_true', help="plot with fixed base width")
    parser.add_argument('--sig_scale', required=False, type=str, default="", help="plot the scaled signal. Supported scalings: [medmad, znorm]")
    # parser.add_argument('--reverse_signal', required=False, action='store_true', help="plot RNA reference/read from 5`-3` and reverse the signal")
    parser.add_argument('--no_pa', required=False, action='store_false', help="skip converting the signal to pA values")
    parser.add_argument('--point_size', required=False, type=int, default=5, help="signal point size [5]")
    parser.add_argument('--base_width', required=False, type=int, default=FIXED_BASE_WIDTH, help="base width when plotting with fixed base width")
    parser.add_argument('--base_shift', required=False, type=int, default=PLOT_BASE_SHIFT, help="the number of bases to shift to align fist signal move")
    parser.add_argument('--plot_limit', required=False, type=int, default=1000, help="limit the number of plots generated")
    parser.add_argument('--sig_plot_limit', required=False, type=int, default=SIG_PLOT_LENGTH, help="maximum number of signal samples to plot")
    parser.add_argument('--stride', required=False, type=int, default=DEFAULT_STRIDE, help="stride used in basecalling network")
    parser.add_argument('--bed', required=False, help="bed file with annotations")
    parser.add_argument('-o', '--output_dir', required=True, help="output dir")
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)


