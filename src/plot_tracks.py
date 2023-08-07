"""
Signal to seQuence alignment Plot - plot
Hiruna Samarakoon - Garvan Medical Institute
hiruna@unsw.edu.au
"""
from bokeh.plotting import output_file, save
from bokeh.layouts import column
import argparse
import os
# import plot
from src import plot_pileup

PLOT_HEIGHT_SCALE = 1000
def run(args):
    print("Info: commands file: {}".format(args.file))
    commands_file = open(args.file, "r")
    num_commands = int(commands_file.readline().strip().split('=')[1])
    print("Info: no. of commands: {}".format(num_commands))

    plot_heights = commands_file.readline().strip().split('=')[1].split(',')
    if plot_heights[0] == '*' or args.auto_height:
        args.auto_height = True
        print("Info: track heights will be set automatically.")
    if not args.auto_height:
        if len(plot_heights) != num_commands:
            raise Exception("Error: number of commands ({}) does not match the number of heights specified ({})".format(num_commands, len(plot_heights)))
        print("Info: specified plot heights: {}".format(str(plot_heights)))

    num_tracks = 0
    pileup = []
    num_plots = []
    total_plots = 0

    for i in range(0, num_commands):
        line_ = commands_file.readline().strip().split()
        tool = str(line_[1])
        # print(tool)
        use_pileup = 0
        if 'plot_pileup' in tool:
            use_pileup = 1
        else:
            raise Exception("Error: {} is not suported in tracks. Please report on github if you want this feature.".format(tool))

        command = line_[2:]
        if '-o' in command:
            idx = command.index("-o")
            del command[idx:idx+2]
        if '--output_dir' in command:
            idx = command.index("--output_dir")
            del command[idx:idx+2]
        if '--return_plot' not in command:
            command.append('--return_plot')

        args_tool = plot_pileup.argparser().parse_args(command)
        p, n_plot = plot_pileup.run(args_tool)
        if p is None or n_plot == 0:
            continue
        p.sizing_mode = 'scale_width'
        num_plots.append(n_plot)
        total_plots += n_plot

        if num_tracks > 0 and args.shared_x:
            if len(p.x_range.js_property_callbacks) > 1:
                raise Exception("Error: more than the expected number of callbacks found ({}). please report on github with a minimal dataset to reproduce this error.".format(len(p.x_range.js_property_callbacks)))
            pileup[num_tracks-1].x_range.js_property_callbacks['change:start'] += p.x_range.js_property_callbacks['change:start']
            p.x_range = pileup[num_tracks-1].x_range
            # p.y_range = pileup[0].y_range
        pileup.append(p)
        # output_file(args.output_dir+"/"+str(num_tracks)+".html", title=num_tracks)
        # save(p)
        num_tracks += 1
    if len(pileup) == 0:
        raise Exception("Error: no plots to plot")
    if args.auto_height:
        j = 0
        for p in pileup:
            # print(num_plots[j],(num_plots[j]/total_plots)*PLOT_HEIGHT_SCALE)
            p.height = int((num_plots[j]/total_plots)*PLOT_HEIGHT_SCALE)
            j += 1
    else:
        j = 0
        for p in pileup:
            p.height = int(plot_heights[i])
            j += 1

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    plot_title = "track"
    if args.tag_name:
        plot_title = args.tag_name
    pileup_output_file_name = args.output_dir + "/" + plot_title + ".html"
    pileup_fig = column(pileup, sizing_mode='stretch_both')
    output_file(pileup_output_file_name, title=plot_title)
    save(pileup_fig)
    print(f'output file: {os.path.abspath(pileup_output_file_name)}')

def argparser():
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument('-f', '--file', required=True, help="commands file")
    parser.add_argument('-o', '--output_dir', required=True, help="output dir")
    parser.add_argument('--tag_name', required=False, type=str, default="", help="a tag name to easily identify the plot")
    parser.add_argument('--shared_x', required=False, action='store_true', help="share x-axis so that all the plots move together")
    parser.add_argument('--auto_height', required=False, action='store_true', help="adjust track height automatically")
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)
