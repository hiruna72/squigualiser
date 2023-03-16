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
import pysam

BASE_LIMIT = 1000
SIG_PLOT_LENGTH = 10000
DEFAULT_STRIDE = 5

BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9


READ_ID = 0
LEN_RAW_SIGNAL = 1
START_RAW = 2
END_RAW = 3
STRAND = 4
READ_ID = 5
LEN_KMER = 6
START_KMER = 7
END_KMER = 8
MATCHES = 9
LEN_KMER = 10
MAPQ = 11

base_color_map = {'A': 'limegreen', 'C': 'blue', 'T': 'red', 'G': 'orange'}

parser = argparse.ArgumentParser()
   
parser.add_argument('-f', '--fasta', required=True, help="fasta file")
parser.add_argument('-r', '--read_id', required=False, help="read id")
parser.add_argument('--base_limit', required=False, type=int, help="maximum number of bases to plot")
parser.add_argument('-s', '--slow5', required=True, help="slow5 file")
parser.add_argument('-a', '--alignment', required=True, help="for read-signal alignment use PAF\nfor reference-signal alignment use SAM/BAM")
parser.add_argument('--region', required=False, type=str, default="", help="[start-end] region to load from SAM/BAM to plot (eg: chr1:6811428-6811467 or chr1:6,811,428-6,811,467)")
parser.add_argument('--tag_name', required=False, type=str, default="", help="a tag name to easily identify the plot")
parser.add_argument('--no_reverse', required=False, action='store_true', help="skip plotting reverse mapped reads")
parser.add_argument('--point_size', required=False, type=int, default=5, help="signal point size [5]")
parser.add_argument('-o', '--output_dir', required=True, help="output dir")

def adjust_before_plotting(ref_seq_len, sam_record, signal_tuple, region_tuple, sig_algn_data, fasta_seq):
    ref_region_start_diff = int(sam_record.reference_start) + 1 - region_tuple[0]
    ref_region_end_diff = int(sam_record.reference_start) + ref_seq_len - region_tuple[1]
    # print("ref_region_start_diff: " + str(ref_region_start_diff))
    if ref_region_end_diff == 0 and ref_region_start_diff == 0:
        return sam_record, signal_tuple, region_tuple, sig_algn_data, fasta_seq

    moves = sig_algn_dic['ss']
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

        sig_algn_dic['ss'] = moves

    return sam_record, signal_tuple, region_tuple, sig_algn_data, fasta_seq

def plot_function(read_id, output_file_name, signal_tuple, sig_algn_data, fasta_tuple, base_limit):
    x = signal_tuple[0]
    x_real = signal_tuple[1]
    y = signal_tuple[2]

    fasta_sequence = fasta_tuple[0]
    query_sequence = fasta_tuple[1]
    cigar_tuple = fasta_tuple[2]

    plot_title = f'{sig_algn_dic["tag_name"]} {read_id}  [signal_start_index,signal_end_index,signal_alignment_start_index:{start_index},{end_index},{sig_algn_data["start_raw"]}]  [seq_length,kmer_start_index,kmer_end_index:{sig_algn_data["len_kmer"]},{sig_algn_data["start_kmer"]},{sig_algn_data["end_kmer"]}]'
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
    base_label_colors = []
    location_plot = 0
    previous_location = location_plot
    initial_location = location_plot
    # draw moves
    moves = sig_algn_data["ss"]

    vlines = []
    base_count = int(sig_algn_data["start_kmer"])
    query_base_count = 0

    num_Is = 0
    num_Ds = 0

    for i in moves:
        previous_location = location_plot
        n_samples = 0
        if 'D' in i:
            i = re.sub('D', '', i)
            n_samples = int(i)
            prev_loc = previous_location
            for j in range(0, n_samples):
                base = fasta_sequence[base_count]
                base_box = BoxAnnotation(left=prev_loc, right=prev_loc+DEFAULT_STRIDE, fill_alpha=0.2, fill_color='white')
                p.add_layout(base_box)

                base_x.append(prev_loc)
                base_y.append(115)
                label = str(base) + "\t" + str(base_count + 1)
                base_label.append(label)
                base_label_colors.append('red')

                prev_loc += DEFAULT_STRIDE
                base_count += 1
                num_Ds += 1

            location_plot = prev_loc

            x = x + list(range(x[-1] + 1, x[-1] + 1 + n_samples * DEFAULT_STRIDE))
            y_add = np.concatenate((y[:previous_location], [0] * n_samples * DEFAULT_STRIDE), axis=0)
            y = np.concatenate((y_add, y[previous_location:]), axis=0)
            x_add = np.concatenate((x_real[:previous_location], [0] * n_samples * DEFAULT_STRIDE), axis=0)
            x_real = np.concatenate((x_add, x_real[previous_location:]), axis=0)

        elif 'I' in i:
            i = re.sub('I', '', i)
            n_samples = int(i)
            num_Is += n_samples
            location_plot += n_samples

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
            base_label_colors.append('black')
            base_count += 1

        if base_count == base_limit:
            break
        if location_plot - initial_location > SIG_PLOT_LENGTH:
            break

    p.renderers.extend(vlines)

    base_annotation = ColumnDataSource(data=dict(base_x=base_x,
                                                 base_y=base_y,
                                                 base_label=base_label,
                                                 colors=base_label_colors))

    base_annotation_labels = LabelSet(x='base_x', y='base_y', text='base_label',
                                      x_offset=5, y_offset=5, source=base_annotation, render_mode='canvas',
                                      text_font_size="9pt", text_color='colors')

    p.add_layout(base_annotation_labels)

    source = ColumnDataSource(data=dict(
        x=x[:location_plot],
        y=y[:location_plot],
        x_real=x_real[:location_plot],
    ))
    p.line('x', 'y', line_width=2, source=source)
    # add a circle renderer with a size, color, and alpha
    p.circle(x[:location_plot], y[:location_plot], size=args.point_size, color="red", alpha=0.5)

    # show the tooltip
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [("x", "@x_real"), ("y", "$y")]
    hover.mode = 'mouse'

    plot_title = f'{sig_algn_dic["tag_name"]} dels:{num_Ds}b ins:{num_Is}samples {read_id}  [signal_start_index,signal_end_index,signal_alignment_start_index:{start_index},{end_index},{sig_algn_data["start_raw"]}]  [seq_length,kmer_start_index,kmer_end_index:{sig_algn_data["len_kmer"]},{sig_algn_data["start_kmer"]},{sig_algn_data["end_kmer"]}]'
    p.title = plot_title

    output_file(output_file_name, title=read_id)
    save(p)

args = parser.parse_args()

if args.fasta:
    print(f'fasta file: {args.fasta}')

if args.base_limit:
    BASE_LIMIT = args.base_limit

print(f'signal file: {args.slow5}')

use_paf = 0
if args.alignment:
    print(f'alignment file: {args.alignment}')
    alignment_extension = args.alignment[-4:]
    if alignment_extension == ".bam" or alignment_extension == ".sam":
        use_paf = 0
    elif alignment_extension == ".paf":
        use_paf = 1
    else:
        print("error please provide the alignment file with correct extension")
        exit()

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

# open signal file
s5 = pyslow5.Open(args.slow5, 'r')

if use_paf == 1:
    with open(args.alignment, "r") as handle:
        fasta_reads = Fasta(args.fasta)
        for paf_record in parse_paf(handle):
            read_id = paf_record.query_name
            fasta_seq = fasta_reads.get_seq(name=read_id, start=1, end=int(paf_record.target_length)).seq
            output_file_name = args.output_dir+"/"+read_id+".html"
            print(f'output file: {output_file_name}')
            print(f'read_id: {read_id}')

            x = []
            x_real = []
            y = []
            base_limit = BASE_LIMIT
            sig_plot_length = SIG_PLOT_LENGTH

            read = s5.get_read(read_id, pA=True, aux=["read_number", "start_mux"])
            if read is not None:
                start_index = paf_record.query_start
                end_index = read['len_raw_signal']

                x = list(range(1, end_index - start_index + 1))
                x_real = list(range(start_index + 1, end_index + 1))  # 1based
                y = read['signal'][start_index:end_index]
            signal_tuple = (x, x_real, y)

            sig_algn_dic = {}
            sig_algn_dic['query_name'] = paf_record.query_name
            sig_algn_dic['start_raw'] = paf_record.query_start
            sig_algn_dic['len_kmer'] = paf_record.target_length
            sig_algn_dic['start_kmer'] = paf_record.target_start
            sig_algn_dic['end_kmer'] = paf_record.target_end
            sig_algn_dic['tag_name'] = args.tag_name

            moves_string = paf_record.tags['ss'][2]
            moves_string = re.sub('D', 'D,', moves_string)
            moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
            moves = re.split(r',+', moves_string)
            sig_algn_dic['ss'] = moves

            fasta_t = (fasta_seq, fasta_seq, fasta_seq)
            plot_function(read_id=read_id, output_file_name=output_file_name, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_tuple=fasta_t, base_limit=base_limit)
else:
    samfile = pysam.AlignmentFile(args.alignment, mode='r')
    fasta_reads = Fasta(args.fasta)
    num_plots = 0
    for sam_record in samfile.fetch():
        if sam_record.is_supplementary:
            continue
        if sam_record.is_reverse and args.no_reverse:
            continue
        # if sam_record.query_name != "76ec70fd-6731-45ab-a820-a3a26e29a035":
        #     continue
        cigar_t = sam_record.cigartuples
        ref_seq_len = 0
        for a in cigar_t:
            cig_op = a[0]
            cig_count = a[1]
            if cig_op == BAM_CMATCH or cig_op == BAM_CDEL or cig_op == BAM_CREF_SKIP or cig_op == BAM_CEQUAL or cig_op == BAM_CDIFF:
                ref_seq_len = ref_seq_len + cig_count

        # print("ref_seq_len: " + str(ref_seq_len))
        if ref_seq_len < BASE_LIMIT:
            base_limit = ref_seq_len
        else:
            base_limit = BASE_LIMIT
        ref_name = ""
        ref_start = 0
        ref_end = 0
        if args.region != "":
            args_region = re.sub(',', '', args.region)
            print(args_region)
            pattern = re.compile("^[a-z]+[0-9]+\:[0-9]+\-[0-9]+")
            if not pattern.match(args_region):
                print("Error: region provided is not in correct format")
                exit(1)
            ref_name = args_region.split(":")[0]
            ref_start = int(args_region.split(":")[1].split("-")[0])
            ref_end = int(args_region.split(":")[1].split("-")[1])

            if ref_name != sam_record.reference_name:
                print("Warning: sam record's reference name [" + sam_record.reference_name + "] and the name specified are different [" + ref_name + "]")
                continue
            # print("ref_start: " + str(ref_start))
            # print("ref_end: " + str(ref_end))
            # print("sam_record.reference_start: " + str(sam_record.reference_start + 1))
            # print("sam_record.reference_start + ref_seq_len: " + str(int(sam_record.reference_start) + ref_seq_len))
            if ref_start >= (int(sam_record.reference_start) + ref_seq_len):
                print("Warning: sam record's region and the region specified do not overlap")
                continue
            elif ref_end <= int(sam_record.reference_start) + 1:
                print("Warning: sam record's region and the region specified do not overlap")
                continue

            if ref_end > (int(sam_record.reference_start) + ref_seq_len):
                ref_end = int(sam_record.reference_start) + ref_seq_len
            if ref_start < (int(sam_record.reference_start) + 1):
                ref_start = int(sam_record.reference_start) + 1

            if (ref_end - ref_start + 1) < BASE_LIMIT:
                base_limit = ref_end - ref_start + 1
        else:
            ref_name = sam_record.reference_name
            ref_start = int(sam_record.reference_start) + 1
            ref_end = int(sam_record.reference_start) + ref_seq_len

        print("plot region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, sam_record.query_name))
        read_id = sam_record.query_name
        fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
        output_file_name = args.output_dir+"/"+read_id+".html"
        print(f'output file: {output_file_name}')

        x = []
        x_real = []
        y = []
        sig_plot_length = SIG_PLOT_LENGTH
        rq_paf = sam_record.get_tag("rq").split(',')

        read = s5.get_read(read_id, pA=True, aux=["read_number", "start_mux"])
        if read is not None:
            # print("read_id:", read['read_id'])
            # print("len_raw_signal:", read['len_raw_signal'])
            start_index = int(rq_paf[START_RAW])
            end_index = read['len_raw_signal']

            x = list(range(1, end_index - start_index + 1))
            x_real = list(range(start_index+1, end_index+1))             # 1based
            y = read['signal'][start_index:end_index]

        if sam_record.is_reverse:
            x_real.reverse()
            y = np.flip(y)

        signal_tuple = (x, x_real, y)
        region_tuple = (ref_start, ref_end)

        sig_algn_dic = {}
        sig_algn_dic['query_name'] = sam_record.query_name
        sig_algn_dic['start_raw'] = rq_paf[START_RAW]
        sig_algn_dic['len_kmer'] = base_limit
        sig_algn_dic['start_kmer'] = rq_paf[START_KMER]
        sig_algn_dic['end_kmer'] = base_limit

        strand_dir = "+"
        if sam_record.is_reverse:
            strand_dir = "-"

        sig_algn_dic['tag_name'] = args.tag_name + " " + ref_name + ":" + str(ref_start) + "-" + str(ref_end) + " (" + strand_dir + ")"

        moves_string = sam_record.get_tag("ss")
        moves_string = re.sub('D', 'D,', moves_string)
        moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
        moves = re.split(r',+', moves_string)

        sig_algn_dic['ss'] = moves

        # print(len(fasta_seq))
        # print(fasta_seq)
        sam_record, signal_tuple, region_tuple, sig_algn_dic, fasta_seq = adjust_before_plotting(ref_seq_len, sam_record, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)

        fasta_t = (fasta_seq, sam_record.query_sequence, sam_record.cigartuples)
        plot_function(read_id=read_id, output_file_name=output_file_name, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_tuple=fasta_t, base_limit=base_limit)
        num_plots += 1

print("Number of plots: {}".format(num_plots))

s5.close()
