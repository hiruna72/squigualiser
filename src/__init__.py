from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from src import sqp, reform
from ._version import __version__

modules = ['sqp', 'reform']
version = "ideal-goggles     {}".format(__version__)

def main():
    parser = ArgumentParser(
        'ideal-goggles',
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