from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, RawTextHelpFormatter
from src import plot, reform, realign, plot_pileup, plot_tracks, calculate_offsets
from ._version import __version__

modules = ['plot', 'reform', 'realign', 'plot_pileup', 'plot_tracks', 'calculate_offsets']
module_help = {
    'plot': {'help': 'Plot read/reference - signal alignments.'},
    'reform': {'help': "Convert basecaller's move table to ss string format."},
    'realign': {'help': 'Realign signal to reference using cigar string and the move table.'},
    'plot_pileup': {'help': 'Plot a reference - signal alignment pileup.'},
    'plot_tracks': {'help': 'Plot multiple reference - signal alignment pileup tracks.'},
    'calculate_offsets': {'help': 'A utility program to calculate the most significant base index given a kmer model or a read - signal alignment.'},
}
version = "squigualiser     {}".format(__version__)

def main():
    parser = ArgumentParser(
        description='squigualiser - A simple tool to Visualise nanopore raw signal-base alignment.',
        epilog='''
See slow5.page.link/squigualiser for a detailed description of these command-line options.

Citation:...
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
    subparsers.required = True
    for module in modules:
        mod = globals()[module]
        p = subparsers.add_parser(module, help=module_help[module]['help'], parents=[mod.argparser()])
        p.set_defaults(func=mod.run)

    args = parser.parse_args()
    args.func(args)
