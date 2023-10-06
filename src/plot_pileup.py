"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
import numpy as np
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import HoverTool, WheelZoomTool, ColumnDataSource, Label, LabelSet, Segment, Arrow, NormalHead, FreehandDrawTool, Range1d, CustomJS
from bokeh.layouts import column
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
import math
from src import bed_annotation
from src import plot_utils
import cProfile, pstats, io
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
# PLOT_HEIGHT = 900
PLOT_Y_MARGIN = 1
SUBPLOT_X = -100
PLOT_BASE_SHIFT = 0
PLOT_X_PADDING = 200

BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK = range(10)
READ_ID, LEN_RAW_SIGNAL, START_RAW, END_RAW, STRAND, SEQUENCE_ID, LEN_KMER, START_KMER, END_KMER, MATCHES, LEN_KMER, MAPQ = range(12)
SI_START_RAW, SI_END_RAW, SI_START_KMER, SI_END_KMER = range(4)

def plot_function_fixed_width_pileup(read_id, signal_tuple, sig_algn_data, fasta_sequence, base_limit, draw_data, p, num_plots, y_shift, y_max, y_min):
    x = signal_tuple[0]
    x_real = signal_tuple[1]
    y = signal_tuple[2]

    label_position = y_min

    base_x = []
    base_y = []
    base_label = []
    base_label_colors = []

    move_x = []
    move_y = []
    move_label = []
    move_label_colors = []

    sample_label_colors_match = []
    sample_label_colors_insert = []
    location_plot = 0

    x_coordinate = 0
    initial_x_coordinate = x_coordinate

    # draw moves
    moves = sig_algn_data["ss"]
    base_index = sig_algn_data["start_kmer"]
    num_Is = 0
    num_Ds = 0
    fixed_width_x = [0.0]
    line_segment_x = []

    base_box_details = {'left': [], 'right': [], 'fill_color': []}
    flag_base_index_bound = 0
    num_samples_in_insertion = 0

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
                if num_plots == -1 or draw_data['overlap_only']:
                    base_label_colors.append('black')
                    base_box_details['left'].append(prev_loc)
                    base_box_details['right'].append(prev_loc + draw_data["fixed_base_width"])
                    base_box_details['fill_color'].append(plot_utils.get_base_color_map()[base])
                else:
                    base_label_colors.append('red')
                    base_box_details['left'].append(prev_loc)
                    base_box_details['right'].append(prev_loc + draw_data["fixed_base_width"])
                    base_box_details['fill_color'].append('white')

                base_x.append(prev_loc)
                base_y.append(label_position+y_shift)
                label = str(base) + "\n" + str(base_index + 1)
                base_label.append(label)

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
                sample_label_colors_match.append('None')
                sample_label_colors_insert.append('None')

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
                sample_label_colors_match.append('None')
                sample_label_colors_insert.append('purple')
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
            base_box_details['fill_color'].append(plot_utils.get_base_color_map()[base])
            line_segment_x.append(location_plot)
            if num_plots == -1 or draw_data['overlap_only']:
                base_x.append(previous_location)
                base_y.append(label_position + y_shift)
                label = str(base) + "\n" + str(base_index + 1)
                base_label.append(label)
                base_label_colors.append('black')
            else:
                if num_samples_in_insertion > 0:
                    move_x.append(previous_location)
                    move_y.append(label_position + y_shift)
                    move_label.append(str(num_samples_in_insertion))
                    move_label_colors.append('purple')
                move_x.append(previous_location+draw_data["fixed_base_width"]/3)
                move_y.append(label_position + y_shift)
                move_label.append(str(n_samples))
                move_label_colors.append('black')

            for j in range(0, n_samples - num_samples_in_insertion):
                sample_label_colors_match.append('red')
                sample_label_colors_insert.append('None')
            base_index += 1
            num_samples_in_insertion = 0

        if base_index - sig_algn_data["start_kmer"] == base_limit:
            break
        if x_coordinate - initial_x_coordinate > draw_data["sig_plot_limit"]:
            break

    line_segment_source = ColumnDataSource(dict(x=line_segment_x, x1=line_segment_x, y=[y_min+y_shift] * len(line_segment_x), y1=[y_max+y_shift] * len(line_segment_x)))
    glyph = Segment(x0="x", y0="y", x1="x1", y1="y1", line_color="#8b4513", line_width=1, line_alpha=0.5)

    base_annotation = ColumnDataSource(data=dict(base_x=base_x, base_y=base_y, base_label=base_label, colors=base_label_colors))
    base_annotation_labels = LabelSet(x='base_x', y='base_y', text='base_label', x_offset=5, y_offset=5, source=base_annotation, text_font_size="7pt", text_color='colors')

    move_annotation = ColumnDataSource(data=dict(move_x=move_x, move_y=move_y, move_label=move_label, colors=move_label_colors))
    move_annotation_labels = LabelSet(x='move_x', y='move_y', text='move_label', x_offset=5, y_offset=5, source=move_annotation, text_font_size="5pt", text_color='colors', visible=True)

    fixed_width_x = fixed_width_x[1:]
    source = ColumnDataSource(data=dict(x=fixed_width_x[:x_coordinate], y=y[:x_coordinate]+y_shift, x_real=x_real[:x_coordinate], y_real=y[:x_coordinate]))

    if (draw_data['overlap_only'] and num_plots == 0) or not draw_data['overlap_only']:
        p.quad(top=y_max+y_shift, bottom=y_min+y_shift, left=base_box_details['left'], right=base_box_details['right'], color=base_box_details['fill_color'], alpha=0.75)
        p.add_glyph(line_segment_source, glyph)
        p.add_layout(base_annotation_labels)
        if draw_data["plot_num_samples"]:
            p.add_layout(move_annotation_labels)
            x_callback_move_annotation = CustomJS(args=dict(move_annotation_labels=move_annotation_labels, init_font_size=move_annotation_labels.text_font_size[:-2], init_xrange=PLOT_X_RANGE), code="""
            let xzoom = (init_font_size * init_xrange) / (cb_obj.end - cb_obj.start);
            move_annotation_labels['text_font_size'] = String(xzoom) + 'pt';
            """)
            p.x_range.js_on_change('start', x_callback_move_annotation)

    if num_plots != -1:
        if num_plots == 0:
            sample_label_colors_insert[0] = 'purple'
            sample_label_colors_match[0] = 'red'
        sample_labels_match = p.circle(fixed_width_x[:x_coordinate], y[:x_coordinate]+y_shift, radius=draw_data["point_size"], color=sample_label_colors_match, alpha=0.5, legend_label='match', visible=False)
        sample_labels_insert = p.circle(fixed_width_x[:x_coordinate], y[:x_coordinate]+y_shift, radius=draw_data["point_size"], color=sample_label_colors_insert, alpha=0.5, legend_label='insertion', visible=False)

    if not draw_data['overlap_only']:
        if num_plots != -1:
            p.line('x', 'y', name="sig_plot_line", line_width=2, source=source)
    if not draw_data['no_overlap']:
        leg_lable = "read_" + str(num_plots+1)
        if num_plots != -1:
            p.line('x', 'y_real', line_width=2, source=source, legend_label=leg_lable)

    # show the tooltip
    hover = p.select(dict(type=HoverTool))
    hover.renderers = p.select(name="sig_plot_line")
    hover.tooltips = [("x", "@x_real"), ("y", "@y_real")]
    hover.mode = 'mouse'

    signal_region = ""
    indels = f'deletions(bases): {num_Ds} insertions(samples): {num_Is}'
    signal_region = f'[{int(x_real[0])}-{int(x_real[x_coordinate - 1])}]'
    y_plot = y[:x_coordinate]+y_shift
    y_median = np.nanmedian(y_plot)
    y_max = np.nanmax(y_plot)
    y_min = np.nanmin(y_plot)
    arrow = Arrow(end=NormalHead(fill_color="orange", size=10), x_start=-2, y_start=y_median, x_end=-1, y_end=y_median)
    subplot_labels = None
    if num_plots != -1 and not draw_data['overlap_only']:
        sub_plot_y_shift = (y_max - y_min)/6
        source_subplot_labels = ColumnDataSource(data=dict(x=[SUBPLOT_X, SUBPLOT_X, SUBPLOT_X, SUBPLOT_X], y=[y_median+sub_plot_y_shift*1, y_median, y_median-sub_plot_y_shift*1, y_median-sub_plot_y_shift*2], tags=[signal_region, indels, read_id, "read: "+str(num_plots+1)]))
        subplot_labels = LabelSet(x='x', y='y', text='tags', text_font_size="7pt", x_offset=5, y_offset=5, source=source_subplot_labels)
        p.add_layout(subplot_labels)
        p.add_layout(arrow)
    elif num_plots == -1 or (draw_data['overlap_only'] and num_plots == 0):
        source_subplot_labels = ColumnDataSource(data=dict(x=[SUBPLOT_X], y=[y_median], tags=["overlap"]))
        subplot_labels = LabelSet(x='x', y='y', text='tags', text_font_size="7pt", x_offset=5, y_offset=5, source=source_subplot_labels)
        p.add_layout(subplot_labels)
        p.add_layout(arrow)
    if subplot_labels:
        x_callback_subplot_labels = CustomJS(args=dict(subplot_labels=subplot_labels, init_font_size=subplot_labels.text_font_size[:-2], init_xrange=PLOT_X_RANGE), code="""
        let xzoom = (init_font_size * init_xrange) / (cb_obj.end - cb_obj.start);
        subplot_labels['text_font_size'] = String(xzoom) + 'pt';
        """)
        p.x_range.js_on_change('start', x_callback_subplot_labels)
    
    x_callback_base_annotation = CustomJS(args=dict(base_annotation_labels=base_annotation_labels, init_font_size=base_annotation_labels.text_font_size[:-2], init_xrange=PLOT_X_RANGE), code="""
    let xzoom = (init_font_size * init_xrange) / (cb_obj.end - cb_obj.start);
    base_annotation_labels['text_font_size'] = String(xzoom) + 'pt';
    """)
    p.x_range.js_on_change('start', x_callback_base_annotation)

    return p, location_plot, base_index
def run(args):
    if args.list_profile:
        plot_utils.list_profiles_base_shift()
        return
    else:
        if args.file == "":
            raise Exception("Error: the following argument is required: -f/--file")
        if args.slow5 == "":
            raise Exception("Error: the following argument is required: -s/--slow5")
        if args.alignment == "":
            raise Exception("Error: the following argument is required: -a/--alignment")
        if args.region == "":
            raise Exception("Error: the following argument is required: --region")
        if not args.return_plot and args.output_dir == "":
            raise Exception("Error: one of the arguments -o/--output_dir --return_plot is required")

    if args.cprofile:
        pr = cProfile.Profile()
        pr.enable()



    use_fasta = 0
    if args.file:
        print(f'sequence file: {args.file}')
        if args.file[-6:] == ".fastq" or args.file[-6:] == ".fq.gz" or args.file[-3:] == ".fq":
            use_fasta = 0
        elif args.file[-6:] == ".fasta" or args.file[-3:] == ".fa":
            use_fasta = 1
        else:
            raise Exception("Error: please provide the sequence file with correct extension")

    max_location_plot = 0
    max_base_index = 0
    plot_sig_ref_flag = 0
    use_paf = 0
    if args.alignment:
        print(f'alignment file: {args.alignment}')
        alignment_extension = args.alignment[-4:]
        if alignment_extension == ".bam" or alignment_extension == ".sam":
            use_paf = 0
            plot_sig_ref_flag = 1
        elif alignment_extension == ".paf":
            raise Exception("Error: For signal to reference alignment the paf file must be bgzip compressed and tabix indexed")
        elif args.alignment[-7:] == ".paf.gz":
            use_paf = 1
            plot_sig_ref_flag = 1
            index_file = args.alignment + ".tbi"
            if not os.path.exists(index_file):
                raise Exception("Error: please provide a bgzip compressed and tabix indexed paf file to extract the regions")
        else:
            raise Exception("Error: please provide the alignment file with correct extension")

    if plot_sig_ref_flag == 0:
        raise Exception("Error: If you are trying to plot just reads please use plot.py. Otherwise provide the correct signal-reference alignment.")

    if use_paf == 0 and use_fasta == 0:
        print("please provide a .fasta or .fa file when using SAM/BAM")

    # if not args.fixed_width:
    #     raise Exception("Error: pileup works only with fixed base width. Provide the argument --fixed_width")

    if args.base_limit:
        base_limit = args.base_limit
    else:
        base_limit = BASE_LIMIT
    print(f'signal file: {args.slow5}')

    if args.read_id != "":
        args.plot_limit = 1
    if args.read_list != "":
        print(f'read_id list file: {args.read_list}')
        read_id_list = list(line.strip() for line in open(args.read_list))

    bed_dic = {}
    if args.bed:
        print(f'bed file: {args.bed}')
        bed_dic = bed_annotation.create_bed_dic(args)
    if args.plot_reverse:
        print("Info: reads mapped to the reverse strand will be plotted")
    if args.return_plot:
        print("Info: plots will be returned without saving")
    else:
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
    draw_data["plot_y_margin"] = PLOT_Y_MARGIN
    draw_data["no_overlap"] = args.no_overlap
    draw_data["overlap_only"] = args.overlap_only
    draw_data["plot_num_samples"] = args.plot_num_samples
    draw_data["bed_labels"] = args.print_bed_labels

    if args.profile == "":
        draw_data["base_shift"] = args.base_shift
    else:
        if args.plot_reverse:
            draw_data["base_shift"] = plot_utils.search_for_profile_base_shift(args.profile)[1]
        else:
            draw_data["base_shift"] = plot_utils.search_for_profile_base_shift(args.profile)[0]

    kmer_correction = 0
    if args.profile != "":
        kmer_correction = -1*(plot_utils.search_for_profile_base_shift(args.profile)[0] + plot_utils.search_for_profile_base_shift(args.profile)[1])

    sig_algn_dic = {}

    y_shift = 0
    prev_y_max = 0
    prev_y_min = 0
    p = plot_utils.create_figure(args, plot_mode=1)
    previous_plot = p

    if use_paf == 1 and plot_sig_ref_flag == 0:
        raise Exception("Error: If you are trying to plot just reads please use plot.py. Otherwise provide the correct signal-reference alignment.")
    elif use_paf == 0 and plot_sig_ref_flag == 1: # using sam/bam
        print("Info: Signal to reference method using SAM/BAM ...")
        fasta_reads = Fasta(args.file)

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
        # if args.bed:
        #     p = bed_annotation.plot_bed_annotation(p=p, ref_id=args_ref_name, bed_dic=bed_dic, sig_algn_data=sig_algn_dic, draw_data=draw_data, base_limit=base_limit)

        for sam_record in samfile.fetch(contig=args_ref_name, start=args_ref_start, stop=args_ref_end):
            if args_ref_name != sam_record.reference_name:
                raise Exception("Error: sam record's reference name [" + sam_record.reference_name + "] and the name specified are different [" + ref_name + "]")
            read_id = sam_record.query_name
            if sam_record.is_supplementary or sam_record.is_unmapped or sam_record.is_secondary:
                continue
            if args.plot_reverse is True and sam_record.is_reverse is False:
                continue
            if args.plot_reverse is False and sam_record.is_reverse is True:
                continue
            if args.read_id != "" and read_id != args.read_id:
                continue
            if args.read_list != "" and read_id not in read_id_list:
                continue
            if not sam_record.has_tag("ss"):
                raise Exception("Error: ss string is missing for the read_id {} in {}".format(read_id, args.alignment))

            ref_seq_len = 0
            data_is_rna = 0
            start_index = -1
            if sam_record.has_tag("si"):
                si_tag = sam_record.get_tag("si").split(',')
                start_index = int(si_tag[SI_START_RAW])
                end_index = int(si_tag[SI_END_RAW])
                ref_seq_len = int(si_tag[SI_END_KMER]) - int(si_tag[SI_START_KMER])
                reference_start = int(si_tag[SI_START_KMER])
                if int(si_tag[SI_START_KMER]) > int(si_tag[SI_END_KMER]):  # if RNA start_kmer>end_kmer in paf
                    data_is_rna = 1
                    if not args.rna:
                        print("Info: data is detected as RNA")
                        raise Exception("Error: data is not specified as RNA. Please provide the argument --rna ")
                    ref_seq_len = int(si_tag[SI_START_KMER]) - int(si_tag[SI_END_KMER])
                    reference_start = int(si_tag[SI_END_KMER]) + kmer_correction

            else:
                raise Exception("Error: sam record does not have a 'si' tag.")
            # print("ref_seq_len: " + str(ref_seq_len))
            if ref_seq_len < BASE_LIMIT:
                base_limit = ref_seq_len
            else:
                base_limit = BASE_LIMIT
            sam_record_reference_end = reference_start + ref_seq_len #1based closed
            if True: # if not args.loose_bound:
                if data_is_rna == 1:
                    if args_ref_start < reference_start + 1 - kmer_correction:
                        continue
                else:
                    if args_ref_start < reference_start + 1:
                        continue
                if args_ref_end > sam_record_reference_end:
                    continue

            ref_name = sam_record.reference_name
            ref_start = reference_start + 1
            ref_end = ref_start + base_limit - 1 #ref_end is 1based closed
            if args_ref_start > ref_start:
                ref_start = args_ref_start
                if (ref_start + base_limit - 1) < sam_record_reference_end:
                    ref_end = ref_start + base_limit - 1
                else:
                    ref_end = sam_record_reference_end
            if args_ref_end < ref_end:
                ref_end = args_ref_end
            if ref_end < ref_start:
                print("Warning: a corner case has hit because the kmer_length used is larger than 1. This alignment will be skipped")
                continue
            base_limit = ref_end - ref_start + 1

            # print("ref_start: {}".format(ref_start))
            # print("ref_end: {}".format(ref_end))
            # print("ref_seq_len: {}".format(ref_seq_len))
            # print("{}".format(sam_record.cigarstring))
            # print("base_limit: {}".format(base_limit))
            if sam_record.is_reverse and data_is_rna == 1:
                raise Exception("Error: A transcript is  always sequenced from 3` to 5`. Squigualiser only supports reads mapped to the transcriptome.")

            if data_is_rna == 1:
                print("plot (RNA 3'->5') region: {}:{}-{}\tread_id: {}".format(ref_name, ref_end, ref_start, read_id))
                fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                fasta_seq = fasta_seq.upper()
            else:
                if sam_record.is_reverse:
                    print("plot (-) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                    fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                    fasta_seq = fasta_seq.upper()
                    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    # fasta_seq = "".join(nn[n] for n in reversed(fasta_seq))
                    fasta_seq = "".join(nn[n] for n in fasta_seq)
                else:
                    print("plot (+) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                    fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                    fasta_seq = fasta_seq.upper()
            if not bool(re.match('^[ACGTUMRWSYKVHDBN]+$', fasta_seq)):
                raise Exception("Error: base characters other than A,C,G,T/U,M,R,W,S,Y,K,V,H,D,B,N were detected. Please check your sequence files.")

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

            scaling_str = "no scaling"
            if args.sig_scale == "medmad" or args.sig_scale == "znorm" or args.sig_scale == "scaledpA":
                scaling_str = args.sig_scale
                draw_data["plot_y_margin"] = 0.1
            elif not args.sig_scale == "":
                raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))

            scale_params = {}
            if args.sig_scale == "scaledpA":
                if not args.no_pa:
                    raise Exception("Error: given --sig_scale method: {} required the signal to be converted to pA levels. Please remove --no_pa argument".format(args.sig_scale))
                for tag in ["sc", "sh"]:
                    if sam_record.has_tag(tag):
                        scale_params[tag] = sam_record.get_tag(tag)
                    else:
                        raise Exception("Error: given --sig_scale method: {} requires {} tag in the alignment file".format(args.sig_scale, tag))
            y = plot_utils.scale_signal(y, args.sig_scale, scale_params)

            moves_string = sam_record.get_tag("ss")
            moves_string = re.sub('D', 'D,', moves_string)
            moves_string = re.sub('I', 'I,', moves_string).rstrip(',')
            moves = re.split(r',+', moves_string)

            if data_is_rna == 0:
                strand_dir = "(DNA +)"
                if sam_record.is_reverse:
                    strand_dir = "(DNA -)"
                    x_real.reverse()
                    y = np.flip(y)
                    moves.reverse()
            if data_is_rna == 1:
                strand_dir = "(RNA 3'->5')"
                fasta_seq = fasta_seq[::-1]
                if sam_record.is_reverse:
                    raise Exception("Error: data is rna and sam record is reverse mapped. This is not implemented yet. Please report")

            signal_tuple = (x, x_real, y)
            region_tuple = (ref_start, ref_end, reference_start, reference_start+ref_seq_len)

            sig_algn_dic['start_kmer'] = 0
            sig_algn_dic['ref_start'] = ref_start
            sig_algn_dic['ref_end'] = ref_end
            sig_algn_dic['pa'] = args.no_pa
            sig_algn_dic['plot_sig_ref_flag'] = plot_sig_ref_flag
            sig_algn_dic['data_is_rna'] = data_is_rna

            sig_algn_dic['ss'] = moves
            # print(len(moves))
            # print(fasta_seq)
            signal_tuple, region_tuple, sig_algn_dic, fasta_seq = plot_utils.adjust_before_plotting(ref_seq_len, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)
            # print(len(sig_algn_dic['ss']))
            # if args.auto_base_shift and num_plots == 0:
            #     draw_data["base_shift"] = plot_utils.calculate_base_shift(signal_tuple[2], fasta_seq, sig_algn_dic['ss'], args)
            #     print("automatically calculated base_shift: {}".format(draw_data["base_shift"]))

            if draw_data["base_shift"] < 0:
                abs_base_shift = abs(draw_data["base_shift"])
                x = signal_tuple[0]
                x_real = signal_tuple[1]
                y = signal_tuple[2]
                y_prefix = [np.nan] * abs_base_shift * draw_data["fixed_base_width"]
                y = np.concatenate((y_prefix, y), axis=0)
                x_real = np.concatenate(([1] * abs_base_shift * draw_data["fixed_base_width"], x_real), axis=0)
                x = list(range(1, len(x) + 1 + abs_base_shift * draw_data["fixed_base_width"]))
                signal_tuple = (x, x_real, y)
                moves_prefix = [str(draw_data["fixed_base_width"])] * abs_base_shift
                sig_algn_dic['ss'] = moves_prefix + sig_algn_dic['ss']

            sig_algn_dic['tag_name'] = args.tag_name + indt + "base_shift: " + str(draw_data["base_shift"]) + indt + "scale:" + scaling_str + indt + "fixed_width: " + str(args.base_width) + indt + strand_dir + indt + "region: " + ref_name + ":"

            y_min = math.floor(np.nanmin(y))
            y_max = math.ceil(np.nanmax(y))

            if num_plots == 0 and args.bed and args.overlap_bottom is False:
                draw_data['y_min'] = y_min
                draw_data['y_max'] = y_max
                draw_data['fixed_width'] = True
                p = bed_annotation.plot_bed_annotation(p=p, ref_id=args_ref_name, bed_dic=bed_dic, sig_algn_data=sig_algn_dic, draw_data=draw_data, base_limit=base_limit)

            if not args.no_overlap and not args.overlap_only:
                if num_plots == 0:
                    p, location_plot, base_index = plot_function_fixed_width_pileup(read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data, p=previous_plot, num_plots=-1, y_shift=y_shift, y_min=y_min, y_max=y_max)
                    previous_plot = p
                    prev_y_max = y_max
                    prev_y_min = y_min
                    if max_base_index < base_index:
                        max_base_index = base_index
                    if max_location_plot < location_plot:
                        max_location_plot = location_plot

            if args.overlap_bottom:
                # y_shift = y_shift + prev_y_min + prev_y_max - prev_y_min + draw_data["plot_y_margin"] - y_min
                y_shift = y_shift + prev_y_max + draw_data["plot_y_margin"] - y_min
            else:
                y_shift = y_shift + prev_y_min - draw_data["plot_y_margin"] - y_max
            if args.overlap_only:
                y_shift = 0
            p, location_plot, base_index = plot_function_fixed_width_pileup(read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data, p=previous_plot, num_plots=num_plots, y_shift=y_shift, y_min=y_min, y_max=y_max)
            previous_plot = p
            prev_y_max = y_max
            prev_y_min = y_min
            if max_base_index < base_index:
                max_base_index = base_index
            if max_location_plot < location_plot:
                max_location_plot = location_plot
            num_plots += 1
            if num_plots == args.plot_limit:
                break
        if args.bed and args.overlap_bottom is True:
            draw_data['y_min'] = prev_y_min + y_shift
            draw_data['y_max'] = prev_y_max + y_shift
            draw_data['fixed_width'] = True
            p = bed_annotation.plot_bed_annotation(p=p, ref_id=args_ref_name, bed_dic=bed_dic, sig_algn_data=sig_algn_dic, draw_data=draw_data, base_limit=base_limit)

    elif use_paf == 1 and plot_sig_ref_flag == 1:
        print("Info: Signal to reference method using PAF ...")
        fasta_reads = Fasta(args.file)
        tbxfile = pysam.TabixFile(args.alignment)

        args_region = re.sub(',', '', args.region)
        # print(args_region)
        # pattern = re.compile("^[a-z]+[0-9]+\:[0-9]+\-[0-9]+")
        pattern = re.compile("^.*\:[0-9]+\-[0-9]+")
        if not pattern.match(args_region):
            raise Exception("Error: region provided is not in correct format")
        args_ref_name = args_region.split(":")[0]
        args_ref_start = int(args_region.split(":")[1].split("-")[0])
        args_ref_end = int(args_region.split(":")[1].split("-")[1])
        for paf_record in tbxfile.fetch(args_ref_name, args_ref_start, args_ref_end, parser=pysam.asTuple()):
            # if paf_record[READ_ID] == paf_record[SEQUENCE_ID]:
            #     raise Exception("Error: this paf file is a signal to read mapping.")
            if args_ref_name != paf_record[SEQUENCE_ID]:
                raise Exception("Error: sam record's reference name [" + paf_record[SEQUENCE_ID] + "] and the name specified are different [" + ref_name + "]")
            read_id = paf_record[READ_ID]
            if args.read_id != "" and read_id != args.read_id:
                continue
            if args.read_list != "" and read_id not in read_id_list:
                continue
            if args.plot_reverse is True and paf_record[STRAND] == "+":
                continue
            if args.plot_reverse is False and paf_record[STRAND] == "-":
                continue
            moves_string = ""
            for i in range(12, len(paf_record)):
                tag = paf_record[i].split(':')[0]
                if tag == "ss":
                    moves_string = paf_record[i].split(':')[2]
            if moves_string == "":
                raise Exception("Error: ss string is missing for the read_id {} in {}".format(read_id, args.alignment))

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
                reference_start = int(paf_record[END_KMER]) + kmer_correction
            # print("ref_seq_len: " + str(ref_seq_len))
            if ref_seq_len < BASE_LIMIT:
                base_limit = ref_seq_len
            else:
                base_limit = BASE_LIMIT
            paf_record_reference_end = reference_start + ref_seq_len #1based closed
            if True: # if not args.loose_bound:
                if data_is_rna == 1:
                    if args_ref_start < reference_start + 1 - kmer_correction:
                        continue
                else:
                    if args_ref_start < reference_start + 1:
                        continue
                if args_ref_end > paf_record_reference_end:
                    continue

            ref_name = paf_record[SEQUENCE_ID]
            ref_start = reference_start + 1
            ref_end = ref_start + base_limit - 1 #ref_end is 1based closed
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
                raise Exception("Error: A transcript is  always sequenced from 3` to 5`. Squigualiser only supports reads mapped to the transcriptome.")
            if data_is_rna == 1:
                print("plot (RNA 3'->5') region: {}:{}-{}\tread_id: {}".format(ref_name, ref_end, ref_start, read_id))
                fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                fasta_seq = fasta_seq.upper()
            else:
                if record_is_reverse:
                    print("plot (-) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                    fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                    fasta_seq = fasta_seq.upper()
                    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    # fasta_seq = "".join(nn[n] for n in reversed(fasta_seq))
                    fasta_seq = "".join(nn[n] for n in fasta_seq)
                else:
                    print("plot (+) region: {}:{}-{}\tread_id: {}".format(ref_name, ref_start, ref_end, read_id))
                    fasta_seq = fasta_reads.get_seq(name=ref_name, start=ref_start, end=ref_end).seq
                    fasta_seq = fasta_seq.upper()
            if not bool(re.match('^[ACGTUMRWSYKVHDBN]+$', fasta_seq)):
                raise Exception("Error: base characters other than A,C,G,T/U,M,R,W,S,Y,K,V,H,D,B,N were detected. Please check your sequence files.")

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

            scaling_str = "no scaling"
            if args.sig_scale == "medmad" or args.sig_scale == "znorm" or args.sig_scale == "scaledpA":
                scaling_str = args.sig_scale
                draw_data["plot_y_margin"] = 0.1
            elif not args.sig_scale == "":
                raise Exception("Error: given --sig_scale method: {} is not supported".format(args.sig_scale))
            scale_params = {}
            if args.sig_scale == "scaledpA":
                if not args.no_pa:
                    raise Exception("Error: given --sig_scale method: {} required the signal to be converted to pA levels. Please remove --no_pa argument".format(args.sig_scale))
                for tag in ["sc", "sh"]:
                    for i in range(12, len(paf_record)):
                        if tag == paf_record[i].split(':')[0]:
                            scale_params[tag] = float(paf_record[i].split(':')[2])
                    if tag not in scale_params:
                        raise Exception("Error: required tag '{}' for given --sig_scale method: {} is not found in the alignment file".format(tag, args.sig_scale))
            y = plot_utils.scale_signal(y, args.sig_scale, scale_params)

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
                strand_dir = "(RNA 3'->5')"
                fasta_seq = fasta_seq[::-1]

            signal_tuple = (x, x_real, y)
            region_tuple = (ref_start, ref_end, reference_start, reference_start + ref_seq_len)

            sig_algn_dic['start_kmer'] = 0
            sig_algn_dic['ref_start'] = ref_start
            sig_algn_dic['ref_end'] = ref_end
            sig_algn_dic['pa'] = args.no_pa
            sig_algn_dic['plot_sig_ref_flag'] = plot_sig_ref_flag
            sig_algn_dic['data_is_rna'] = data_is_rna

            sig_algn_dic['ss'] = moves
            # print(len(moves))
            # print(fasta_seq)
            signal_tuple, region_tuple, sig_algn_dic, fasta_seq = plot_utils.adjust_before_plotting(ref_seq_len, signal_tuple, region_tuple, sig_algn_dic, fasta_seq)

            # if args.auto_base_shift and num_plots == 0:
            #     draw_data["base_shift"] = plot_utils.calculate_base_shift(signal_tuple[2], fasta_seq, sig_algn_dic['ss'], args)
            #     print("automatically calculated base_shift: {}".format(draw_data["base_shift"]))

            if draw_data["base_shift"] < 0:
                abs_base_shift = abs(draw_data["base_shift"])
                x = signal_tuple[0]
                x_real = signal_tuple[1]
                y = signal_tuple[2]
                y_prefix = [np.nan] * abs_base_shift * draw_data["fixed_base_width"]
                y = np.concatenate((y_prefix, y), axis=0)
                x_real = np.concatenate(([1] * abs_base_shift * draw_data["fixed_base_width"], x_real), axis=0)
                x = list(range(1, len(x) + 1 + abs_base_shift * draw_data["fixed_base_width"]))
                signal_tuple = (x, x_real, y)
                moves_prefix = [str(draw_data["fixed_base_width"])] * abs_base_shift
                sig_algn_dic['ss'] = moves_prefix + sig_algn_dic['ss']

            sig_algn_dic['tag_name'] = args.tag_name + indt + "base_shift: " + str(draw_data["base_shift"]) + indt + "scale:" + scaling_str + indt + "fixed_width: " + str(args.base_width) + indt + strand_dir + indt + "region: " + ref_name + ":"

            # print(len(sig_algn_dic['ss']))
            y_min = math.floor(np.nanmin(y))
            y_max = math.ceil(np.nanmax(y))

            if num_plots == 0 and args.bed and args.overlap_bottom is False:
                draw_data['y_min'] = y_min
                draw_data['y_max'] = y_max
                draw_data['fixed_width'] = True
                p = bed_annotation.plot_bed_annotation(p=p, ref_id=args_ref_name, bed_dic=bed_dic, sig_algn_data=sig_algn_dic, draw_data=draw_data, base_limit=base_limit)

            if not args.no_overlap and not args.overlap_only:
                if num_plots == 0:
                    p, location_plot, base_index = plot_function_fixed_width_pileup(read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data, p=previous_plot, num_plots=-1, y_shift=y_shift, y_min=y_min, y_max=y_max)
                    previous_plot = p
                    prev_y_max = y_max
                    prev_y_min = y_min
                    if max_base_index < base_index:
                        max_base_index = base_index
                    if max_location_plot < location_plot:
                        max_location_plot = location_plot

            if args.overlap_bottom:
                # y_shift = y_shift + prev_y_min + prev_y_max - prev_y_min + draw_data["plot_y_margin"] - y_min
                y_shift = y_shift + prev_y_max + draw_data["plot_y_margin"] - y_min
            else:
                y_shift = y_shift + prev_y_min - draw_data["plot_y_margin"] - y_max
            if args.overlap_only:
                y_shift = 0
            p, location_plot, base_index = plot_function_fixed_width_pileup(read_id=read_id, signal_tuple=signal_tuple, sig_algn_data=sig_algn_dic, fasta_sequence=fasta_seq, base_limit=base_limit, draw_data=draw_data, p=previous_plot, num_plots=num_plots, y_shift=y_shift, y_min=y_min, y_max=y_max)
            previous_plot = p
            prev_y_max = y_max
            prev_y_min = y_min
            if max_base_index < base_index:
                max_base_index = base_index
            if max_location_plot < location_plot:
                max_location_plot = location_plot

            num_plots += 1
            if num_plots == args.plot_limit:
                break
        if args.bed and args.overlap_bottom is True:
            draw_data['y_min'] = prev_y_min + y_shift
            draw_data['y_max'] = prev_y_max + y_shift
            draw_data['fixed_width'] = True
            p = bed_annotation.plot_bed_annotation(p=p, ref_id=args_ref_name, bed_dic=bed_dic, sig_algn_data=sig_algn_dic, draw_data=draw_data, base_limit=base_limit)

    else:
        raise Exception("Error: You should not have ended up here. Please check your arguments")

    print("Number of plots: {}".format(num_plots))
    if num_plots == 0:
        print("Squigualiser only plots reads that span across the specified region entirely. Reduce the region interval and double check with IGV : {}".format(num_plots))

    if num_plots > 0:
        if sig_algn_dic["data_is_rna"] == 1:
            sig_dir = " ->"
            plot_title = f'{sig_algn_dic["tag_name"]}[{sig_algn_dic["ref_end"]:,}-{sig_algn_dic["ref_end"] - max_base_index + 1:,}]{indt}num reads:{num_plots}{indt}signal dir:{sig_dir}'
        else:
            if args.plot_reverse:
                sig_dir = " <-"
            else:
                sig_dir = " ->"
            plot_title = f'{sig_algn_dic["tag_name"]}[{sig_algn_dic["ref_start"]:,}-{sig_algn_dic["ref_start"] + max_base_index - 1:,}]{indt}num reads:{num_plots}{indt}signal dir:{sig_dir}'
        p.title = plot_title
        p.legend.click_policy = "hide"
        # p.legend.location = 'top_left'
        p.legend.label_text_font_size = '7pt'
        p.legend.background_fill_alpha = 0.5

        if max_location_plot > y_shift:
            if max_location_plot > PLOT_X_RANGE:
                new_x_range = Range1d(0, PLOT_X_RANGE, bounds=(-1*PLOT_X_PADDING, max_location_plot+PLOT_X_PADDING))
            else:
                new_x_range = Range1d(0, max_location_plot, bounds=(-1*PLOT_X_PADDING, max_location_plot+PLOT_X_PADDING))
            new_x_range.js_property_callbacks = p.x_range.js_property_callbacks
            p.x_range = new_x_range

        renderer = p.multi_line([[1, 1]], [[1, 1]], line_width=4, alpha=0.4, color='black')
        draw_tool = FreehandDrawTool(renderers=[renderer], num_objects=50)
        p.add_tools(draw_tool)
        if args.return_plot:
            return p, num_plots
        else:
            pileup_output_file_name = args.output_dir + "/" + "pileup_" + args.tag_name + ".html"
            output_file(pileup_output_file_name, title="pileup_" + args.tag_name)
            print(f'output file: {os.path.abspath(pileup_output_file_name)}')
            save(p)
    elif num_plots == 0 and args.return_plot:
        return None, num_plots

    s5.close()

    if args.cprofile:
        pr.disable()
        s = io.StringIO()
        sort_by = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sort_by)
        ps.print_stats()
        with open("cprofile_run.log", 'w') as f:
            print(s.getvalue(), file=f)
def argparser():
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    parser.add_argument('-f', '--file', required=False, type=str, default="", help="fasta/fa/fastq/fq/fq.gz sequence file")
    parser.add_argument('-a', '--alignment', required=False, type=str, default="", help="for read-signal alignment use PAF\nfor reference-signal alignment use SAM/BAM")
    parser.add_argument('-s', '--slow5', required=False, type=str, default="", help="slow5 file")
    parser.add_argument('--region', required=False, type=str, default="", help="[start-end] 1-based closed interval region to plot. For SAM/BAM eg: chr1:6811428-6811467 or chr1:6,811,428-6,811,467. For PAF eg:100-200.")
    parser.add_argument('--tag_name', required=False, type=str, default="", help="a tag name to easily identify the plot")
    parser.add_argument('-r', '--read_id', required=False, type=str, default="", help="plot the read with read_id")
    parser.add_argument('-l', '--read_list', required=False, type=str, default="", help="a file with read_ids to plot")
    parser.add_argument('--base_limit', required=False, type=int, help="maximum number of bases to plot")
    parser.add_argument('--plot_reverse', required=False, action='store_true', help="plot only reverse mapped reads")
    parser.add_argument('--plot_num_samples', required=False, action='store_true', help="annotate the number of samples for each move")
    parser.add_argument('--rna', required=False, action='store_true', help="specify for RNA reads")
    parser.add_argument('--sig_scale', required=False, type=str, default="", help="plot the scaled signal. Supported scalings: [medmad, znorm, scaledpA]")
    parser.add_argument('--no_pa', required=False, action='store_false', help="skip converting the signal to pA values")
    parser.add_argument('--overlap_bottom', required=False, action='store_true', help="plot the overlap at the bottom")
    parser.add_argument('--no_overlap', required=False, action='store_true', help="skip plotting the overlap")
    parser.add_argument('--overlap_only', required=False, action='store_true', help="plot only the overlap")
    # parser.add_argument('--loose_bound', required=False, action='store_true', help="also plot alignments not completely within the specified region")
    parser.add_argument('--point_size', required=False, type=int, default=0.5, help="signal point radius [0.5]")
    parser.add_argument('--base_width', required=False, type=int, default=FIXED_BASE_WIDTH, help="base width when plotting with fixed base width")
    parser.add_argument('--base_shift', required=False, type=int, default=PLOT_BASE_SHIFT, help="the number of bases to shift to align fist signal move")
    parser.add_argument('--profile', required=False, default="", type=str, help="determine base_shift using preset values")
    parser.add_argument('--list_profile', action='store_true', help="list the available profiles")
    parser.add_argument('--plot_limit', required=False, type=int, default=1000, help="limit the number of plots generated")
    parser.add_argument('--sig_plot_limit', required=False, type=int, default=SIG_PLOT_LENGTH, help="maximum number of signal samples to plot")
    parser.add_argument('--bed', required=False, help="bed file with annotations")
    parser.add_argument('--print_bed_labels',  required=False, action='store_true', help="draw bed annotations with labels")
    parser.add_argument('--cprofile', required=False, action='store_true', help="run cProfile for benchmarking")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-o', '--output_dir', type=str, default="", help="output dir")
    group.add_argument('--return_plot', action='store_true', help="return plot object without saving to output")
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)

