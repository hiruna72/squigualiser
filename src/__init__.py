from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from src import plot, reform, realign, plot_pileup
from ._version import __version__

modules = ['plot', 'reform', 'realign', 'plot_pileup']
version = "squigualiser     {}".format(__version__)

def main():
    parser = ArgumentParser(
        'squigualiser',
        formatter_class=ArgumentDefaultsHelpFormatter
    )

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
        p = subparsers.add_parser(module, parents=[mod.argparser()])
        p.set_defaults(func=mod.run)

    args = parser.parse_args()
    args.func(args)
