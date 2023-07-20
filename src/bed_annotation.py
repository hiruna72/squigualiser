"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
import numpy as np
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import ColumnDataSource, LabelSet, CustomJS
from bokeh.colors import RGB
import argparse
import re

# ref_start is always 1based closed
# ref_end is always 1based closed
# start_kmer is always 0based closed
# end_kmer is always 0based open

PLOT_X_RANGE = 750

DEFAULT_NUM_BED_COLS = 3
DEFAULT_BED_ANNOTATION_COLOR = (75, 126, 246)
BED_CHROM, BED_CHROM_START, BED_CHROM_END, BED_NAME, BED_SCORE, BED_STRAND, BED_THICK_START, BED_THICK_END, BED_ITEM_RGB, BED_BLOCK_COUNT, BED_BLOCK_SIZES, BLOCK_STARTS = range(12)


def adjust_bed_for_rna(bed_content, sig_algn_data):
    ref_end = sig_algn_data["ref_end"]
    new_bed_content = []
    for i in bed_content:
        bed_chrom_start = int(i[BED_CHROM_START])
        bed_chrom_end = int(i[BED_CHROM_END])
        i[BED_CHROM_START] = ref_end - bed_chrom_end
        i[BED_CHROM_END] = ref_end - bed_chrom_start
        if i[BED_CHROM_END] > 0:
            new_bed_content.append(i)
    return new_bed_content
def draw_bed_annotation(p, bed_content, sig_algn_data, draw_data, base_limit, track_shift, track_height):
    moves = sig_algn_data["ss"]
    base_index = sig_algn_data["start_kmer"]
    ref_start = sig_algn_data["ref_start"]
    if sig_algn_data["data_is_rna"]:
        ref_start = 1
    ref_end = sig_algn_data["ref_end"]

    if sig_algn_data["data_is_rna"]:
        bed_content = adjust_bed_for_rna(bed_content, sig_algn_data)
    else:
        new_bed_content = []
        for i in bed_content:
            if int(i[BED_CHROM_END]) > base_index+ref_start:
                new_bed_content.append(i)
        bed_content = new_bed_content
    if len(bed_content) == 0:
        return p
    bed_content.sort(key=lambda x: x[BED_CHROM_START])
    annotation_box_details = {'left': [], 'right': [], 'fill_color': []}
    annotation_label = []
    annotation_label_x = []
    location_plot = 0
    x_coordinate = 0
    initial_x_coordinate = x_coordinate
    flag_base_index_bound = 0

    bed_index = 0
    bed_region_start = int(bed_content[bed_index][BED_CHROM_START])
    bed_region_end = int(bed_content[bed_index][BED_CHROM_END])

    flag_region_active = 1
    flag_interval_open = 0
    if len(bed_content[bed_index]) >= BED_ITEM_RGB + 1:
        bed_region_color = eval(bed_content[bed_index][BED_ITEM_RGB])
    else:
        bed_region_color = DEFAULT_BED_ANNOTATION_COLOR
    for i in moves:
        previous_location = location_plot
        previous_x_coordinate = x_coordinate
        if 'D' in i:
            i = re.sub('D', '', i)
            n_samples = int(i)
            prev_loc = previous_location
            prev_x_cord = previous_x_coordinate
            for j in range(0, n_samples):
                prev_loc += draw_data["fixed_base_width"]
                prev_x_cord += draw_data["fixed_base_width"]
                #start add annotation
                if flag_interval_open == 0 and flag_region_active == 1 and bed_region_start+1 == base_index+ref_start:
                    annotation_box_details['left'].append(prev_loc-draw_data["fixed_base_width"])
                    flag_interval_open = 1
                if flag_interval_open == 1 and flag_region_active == 1 and bed_region_end == base_index+ref_start:
                    annotation_box_details['right'].append(prev_loc)
                    annotation_box_details['fill_color'].append(bed_region_color)
                    bed_name = "{}-{} ".format(bed_region_start+1, bed_region_end)
                    if sig_algn_data["data_is_rna"]:
                        bed_name = "{}-{} ".format(ref_end-bed_region_start, ref_end-bed_region_end+1)
                    if len(bed_content[bed_index]) >= BED_NAME + 1:
                        bed_name += bed_content[bed_index][BED_NAME]
                    annotation_label.append(bed_name)
                    annotation_label_x.append((annotation_box_details['left'][bed_index]))
                    flag_interval_open = 0
                    flag_region_active = 0

                if flag_region_active == 0 and bed_region_end == base_index+ref_start and bed_index < len(bed_content)-1:
                    bed_index += 1
                    prev_bed_end = bed_region_end
                    bed_region_start = int(bed_content[bed_index][BED_CHROM_START])
                    bed_region_end = int(bed_content[bed_index][BED_CHROM_END])
                    if prev_bed_end > bed_region_start:
                        raise Exception("Error: overlapping regions found in bed file")
                    if len(bed_content[bed_index]) >= BED_ITEM_RGB + 1:
                        bed_region_color = eval(bed_content[bed_index][BED_ITEM_RGB])
                    else:
                        bed_region_color = DEFAULT_BED_ANNOTATION_COLOR
                    flag_region_active = 1
                #end add annotation
                base_index += 1
                if base_index - sig_algn_data["start_kmer"] == base_limit:
                    flag_base_index_bound = 1
                    break
            if flag_base_index_bound == 1:
                break
            location_plot = prev_loc
            x_coordinate = prev_x_cord

        elif 'I' in i:
            i = re.sub('I', '', i)
            n_samples = int(i)
            location_plot += n_samples
            x_coordinate += n_samples
        else:
            n_samples = int(i)
            if draw_data['fixed_width']:
                location_plot += draw_data["fixed_base_width"]
            else:
                location_plot += n_samples
            x_coordinate += n_samples
            #start add annotation
            if flag_interval_open == 0 and flag_region_active == 1 and bed_region_start+1 <= base_index+ref_start:
                annotation_box_details['left'].append(previous_location)
                flag_interval_open = 1
            if flag_interval_open == 1 and flag_region_active == 1 and bed_region_end == base_index+ref_start:
                annotation_box_details['right'].append(location_plot)
                annotation_box_details['fill_color'].append(bed_region_color)
                bed_name = "{}-{} ".format(bed_region_start+1, bed_region_end)
                if sig_algn_data["data_is_rna"]:
                    bed_name = "{}-{} ".format(ref_end-bed_region_start, ref_end-bed_region_end+1)
                if len(bed_content[bed_index]) >= BED_NAME + 1:
                    bed_name += bed_content[bed_index][BED_NAME]
                annotation_label.append(bed_name)
                annotation_label_x.append((annotation_box_details['left'][bed_index]))
                flag_interval_open = 0
                flag_region_active = 0

            if flag_region_active == 0 and bed_region_end == base_index+ref_start and bed_index < len(bed_content)-1:
                bed_index += 1
                prev_bed_end = bed_region_end
                bed_region_start = int(bed_content[bed_index][BED_CHROM_START])
                bed_region_end = int(bed_content[bed_index][BED_CHROM_END])
                if prev_bed_end > bed_region_start:
                    raise Exception("Error: overlapping regions found in bed file")
                if len(bed_content[bed_index]) >= BED_ITEM_RGB + 1:
                    bed_region_color = eval(bed_content[bed_index][BED_ITEM_RGB])
                else:
                    bed_region_color = DEFAULT_BED_ANNOTATION_COLOR
                flag_region_active = 1
            #end add annotation
            base_index += 1

        if base_index - sig_algn_data["start_kmer"] == base_limit:
            break
        if x_coordinate - initial_x_coordinate > draw_data["sig_plot_limit"]:
            break

    if flag_interval_open == 1 and flag_region_active == 1:
        if len(annotation_box_details['left']) != len(annotation_box_details['right']) + 1:
            raise Exception("Error: bed annotation interval dimensions are incorrect. Please report")
        annotation_box_details['right'].append(location_plot)
        annotation_box_details['fill_color'].append(bed_region_color)
        bed_name = "{}- ".format(bed_region_start+1)
        if sig_algn_data["data_is_rna"]:
            bed_name = "{}- ".format(ref_end-bed_region_start)
        if len(bed_content[bed_index]) >= BED_NAME + 1:
            bed_name += bed_content[bed_index][BED_NAME]
        annotation_label.append(bed_name)
        annotation_label_x.append((annotation_box_details['left'][bed_index]))
    # print(annotation_box_details)
    p.quad(top=draw_data['y_max']+track_shift+track_height, bottom=draw_data['y_max']+track_shift, left=annotation_box_details['left'], right=annotation_box_details['right'], color=annotation_box_details['fill_color'], alpha=0.75)

    # print(annotation_label)
    if draw_data["bed_labels"]:
        bed_annotation = ColumnDataSource(data=dict(base_x=annotation_label_x, base_label=annotation_label))
        bed_annotation_labels = LabelSet(x='base_x', y=draw_data['y_max']+track_shift, text='base_label', source=bed_annotation, text_font_size="7pt")
        p.add_layout(bed_annotation_labels)

        x_callback = CustomJS(args=dict(bed_annotation_labels=bed_annotation_labels, init_font_size=bed_annotation_labels.text_font_size[:-2], init_xrange=PLOT_X_RANGE), code="""
        let xzoom = (init_font_size * init_xrange) / (cb_obj.end - cb_obj.start);
        bed_annotation_labels['text_font_size'] = String(xzoom) + 'pt';
        """)
        p.x_range.js_on_change('start', x_callback)

    return p
def plot_bed_annotation(p, ref_id, bed_dic, sig_algn_data, draw_data, base_limit):
    if ref_id in bed_dic:
        plot_height = draw_data['y_max'] - draw_data['y_min']
        track_height = plot_height/5
        # print("track_height: {}", track_height)
        track_count = 0
        for key in bed_dic[ref_id]:
            track_count += 1
            draw_bed_annotation(p, bed_dic[ref_id][key], sig_algn_data, draw_data, base_limit, track_height*track_count, track_height)
    return p
def create_bed_dic(args):
    bed_content = []
    bed_num_cols = DEFAULT_NUM_BED_COLS
    count = 0
    with open(args.bed)as f:
        for line in f:
            bed_list = line.strip().split()
            if len(bed_list) < DEFAULT_NUM_BED_COLS:
                raise Exception("Error: minimum {} columns required in a bed file", bed_num_cols)
            if count == 0:
                bed_num_cols = len(bed_list)
            if len(bed_list) != bed_num_cols:
                raise Exception("Error: Number of columns vary for different bed records")
            count += 1
            bed_content.append(bed_list)
    # print(bed_content)
    bed_dic = {}
    bed_chrom_set = set()
    for i in range(0, len(bed_content)):
        bed_chrom_set.add(bed_content[i][BED_CHROM])
    for i in bed_chrom_set:
        bed_dic[i] = {}
    if len(bed_content) > 0:
        if bed_num_cols >= BED_NAME + 1:
            for i in range(0, len(bed_content)):
                if bed_content[i][BED_NAME] not in bed_dic[bed_content[i][BED_CHROM]]:
                    bed_dic[bed_content[i][BED_CHROM]][bed_content[i][BED_NAME]] = []
                bed_dic[bed_content[i][BED_CHROM]][bed_content[i][BED_NAME]].append(bed_content[i])
        else:
            for i in range(0, len(bed_content)):
                if "DEFAULT" not in bed_dic[bed_content[i][BED_CHROM]]:
                    bed_dic[bed_content[i][BED_CHROM]]["DEFAULT"] = []
                bed_dic[bed_content[i][BED_CHROM]]["DEFAULT"].append(bed_content[i])
    return bed_dic


