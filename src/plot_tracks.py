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

def run(args):
    print(args)
    commands_file = open(args.file, "r")
    num_commands = int(commands_file.readline().strip().split('=')[1])
    print("Info: no. of commands: {}".format(num_commands))
    plot_heights = commands_file.readline().strip().split('=')[1].split(',')
    print("Info: specified plot heights: {}".format(str(plot_heights)))
    num_plots = 0
    # Strips the newline character
    pileup = []
    for i in range(0, num_commands):
        line_ = commands_file.readline().strip().split(' ')
        tool = str(line_[1])
        print(tool)
        use_pileup = 0
        if tool == 'plot_pileup' or tool == 'plot_pileup.py':
            use_pileup = 1
        else:
            print("Error: {} is not suported in tracks. Please report on github if you want this feature.".format(tool))
            exit(1)
        command = line_[2:]
        if '--return_plot' not in command:
            print("Error: please specify --return_plot for each command")
            exit(1)
        # print(command)

        args_tool = plot_pileup.argparser().parse_args(command)
        p = plot_pileup.run(args_tool)
        if p is None:
            continue
        p.height = int(plot_heights[i])
        p.sizing_mode = 'scale_width'

        if num_plots > 0 and args.shared_x:
            p.x_range = pileup[0].x_range
            # p.y_range = pileup[0].y_range
        pileup.append(p)
        # output_file(args.output_dir+"/"+str(num_plots)+".html", title=num_plots)
        # save(p)
        num_plots += 1
    if len(pileup) == 0:
        print("Error: no plots to plot")
        exit(1)
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
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(str(e))
        exit(1)
