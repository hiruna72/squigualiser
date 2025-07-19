import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, RawTextHelpFormatter
from src import plot, reform, realign, plot_pileup, plot_tracks, calculate_offsets, metric, plot_signal
from ._version import __version__

class MyParser(ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help(sys.stderr)
        sys.exit(2)

modules = ['plot', 'reform', 'realign', 'plot_pileup', 'plot_tracks', 'calculate_offsets', 'metric', 'plot_signal']
module_help = {
    'plot': {'help': 'Plot read/reference - signal alignments.'},
    'reform': {'help': "Convert basecaller's move table to ss string format."},
    'realign': {'help': 'Realign signal to reference using cigar string and the move table.'},
    'plot_pileup': {'help': 'Plot a reference - signal alignment pileup.'},
    'plot_tracks': {'help': 'Plot multiple reference - signal alignment pileup tracks.'},
    'calculate_offsets': {'help': 'A utility program to calculate the most significant base index given a kmer model or a read - signal alignment.'},
    'metric': {'help': 'A utility program to calculate some statistics of the signal alignment.'},
    'plot_signal': {'help': 'A utility program to plot signal from a slow5 file.'}
}
version = "squigualiser     {}".format(__version__)

def main():
    parser = MyParser(
        description='squigualiser - A simple tool to Visualise nanopore raw signal-base alignment.',
        epilog='''
See https://slow5.bioinf.science/squigualiser for a detailed description of these command-line options.
    
Citation: Samarakoon, H., Liyanage, K., Ferguson, J.M., Parameswaran, S., Gamaarachchi, H. and Deveson, I.W., 2024. Interactive visualization of nanopore sequencing signal data with Squigualiser. Bioinformatics, 40(8), p.btae501.
               ''',
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        '-v', '--version', action='version',
        version='%(prog)s {}'.format(__version__)
    )

    subparsers = parser.add_subparsers(
        title='subcommands', description='valid commands',
        help='additional help', dest='command'
    )
    # subparsers.required = True
    for module in modules:
        mod = globals()[module]
        p = subparsers.add_parser(module, help=module_help[module]['help'], description=module_help[module]['help'], parents=[mod.argparser()])
        p.set_defaults(func=mod.run)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)


    args.func(args)
