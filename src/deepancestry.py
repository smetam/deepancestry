import argparse

from loader import Loader


VERSION = '0.1'


def load(chromosome, start, end, dir):
    loader = Loader(chromosome, start, end, dir)
    loader.fetch()


def main(args):
    print('well done!')
    pass


if __name__ == '__main__':

    """
    Options:
    - fetch
    - prepare
    - train ? should contain trained model based on 1000genomes 
    - predict
    """
    parser = argparse.ArgumentParser(description='Calculate best markers and their informativeness for specified '
                                                 'segment of the chromosome')

    parser.add_argument('-v', '--version', action='version', version='pylai ' + VERSION)
    subparsers = parser.add_subparsers(
        help="You may run on fetch, prepare or predict mode.", dest='mode',
        description='There are 3 modes in pylai.\n'
                    'For detailed usage of each mode:'
                    '    pylai.py mode -h'
                    '-------------------------------------------------------'
    )

    parser_f = subparsers.add_parser('fetch', help="Load the required chromosome from 1000 genomes project")
    parser_f.add_argument('-c', '--chromosome', action='store', type=str, default='18', help='Chromosome number or X')
    parser_f.add_argument('-s', '--start', action='store', type=int, default=None, help='Start position of chromosome')
    parser_f.add_argument('-e', '--end', action='store', type=int, default=None, help='End position of chromosome')
    parser_f.add_argument('--dir', action='store', type=str, default=None, help='Directory to store chromosome')
    # TODO: add option to load from custom ftp
    # parser_f.add_argument('--source', action='store', type=str, default=None, help='Specify custom ftp of chromosome')

    parser_p = subparsers.add_parser('prepare', help="Prepare data for the chromosome")
    parser_p.add_argument('-n', action='store', type=int, default=10, help='number of markers')
    parser_p.add_argument('-c', '--chromosome', action='store', type=int, default=1, help='Chromosome number')
    parser_p.add_argument('-s', '--start', action='store', type=int, default=1, help='Starting position')
    parser_p.add_argument('-e', '--end', action='store', type=int, default=1000000, help='Ending position')

    parser_i = subparsers.add_parser('predict', help="Predict admixture")
    parser_i.add_argument('-n', action='store', type=int, default=10, help='number of markers')
    parser_i.add_argument('-c', '--chromosome', action='store', type=int, default=1, help='Chromosome number')
    parser_i.add_argument('-s', '--start', action='store', type=int, default=1, help='Starting position')
    parser_i.add_argument('-e', '--end', action='store', type=int, default=1000000, help='Ending position')
    parser_i.add_argument('populations', type=str, nargs='*', default=('GBR', 'FIN'), help='Populations')

    args = parser.parse_args()
    print(args.mode)
    mode = args.mode

    if mode == 'fetch':
        load(args.chromosome, args.start, args.end, args.dir)
    elif mode == 'prepare':
        pass

    main(args)
