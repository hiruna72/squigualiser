from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, RawTextHelpFormatter
from src import plot, reform, realign, plot_pileup, plot_tracks, calculate_offsets
from ._version import __version__

modules = ['plot', 'reform', 'realign', 'plot_pileup', 'plot_tracks', 'calculate_offsets']
module_help = {
    'plot': {'help': 'stuff'},
    'reform': {'help': 'stuff'},
    'realign': {'help': 'stuff'},
    'plot_pileup': {'help': 'stuff'},
    'plot_tracks': {'help': 'stuff'},
    'calculate_offsets': {'help': 'stuff'},
}
version = "squigualiser     {}".format(__version__)

def main():
    parser = ArgumentParser(
        description='squigualiser - A simple tool to Visualise nanopore raw signal-base alignment.',
        epilog='''
See https://github.com/hiruna72/squigualiser/blob/main/docs/commands.md for detailed description of these command-line options.

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
