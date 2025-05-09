import sys
import os
import glob
import argparse

from rich_argparse import ArgumentDefaultsRichHelpFormatter
from melon.utils import logger, get_filename
from melon import GenomeProfiler
from . import __version__, AntibioticResistanceGeneProfiler

## customize formatter
ArgumentDefaultsRichHelpFormatter.styles['argparse.prog'] = 'default'
ArgumentDefaultsRichHelpFormatter.styles['argparse.default'] = 'grey50'
ArgumentDefaultsRichHelpFormatter.styles['argparse.metavar'] = 'grey50'
ArgumentDefaultsRichHelpFormatter.styles['argparse.groups'] = '#E34234'
ArgumentDefaultsRichHelpFormatter.styles['argparse.args'] = 'default'


def cli(argv=sys.argv):
    '''
    Entry point for command line interface.
    '''
    parser = argparse.ArgumentParser(
        description=f'Argo v{__version__}: species-resolved profiling of antibiotic resistance genes in complex metagenomes through long-read overlapping',
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False)

    parser.add_argument(
        dest='files',
        nargs='+',
        metavar='file',
        help='Input fasta <*.fa|*.fasta> or fastq <*.fq|*.fastq> file, gzip optional <*.gz>.')

    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-d',
        '--db',
        metavar='DIR',
        required=True,
        help='Unzipped database folder, should contains <prot.fa|sarg.fa>, <nucl.*.fa|sarg.*.fa> and metadata files.')

    required.add_argument(
        '-o',
        '--output',
        metavar='DIR',
        required=True,
        help='Output folder.')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        '-t',
        '--threads',
        metavar='INT',
        type=int,
        default=os.cpu_count(),
        help=f'Number of threads. [{os.cpu_count()}]')

    optional.add_argument(
        '--plasmid',
        action='store_true',
        help='List ARGs carried by plasmids.')

    optional.add_argument(
        '--skip-melon',
        action='store_true',
        help='Skip Melon for genome copy estimation.')

    optional.add_argument(
        '--skip-clean',
        action='store_true',
        help='Skip cleaning, keep all temporary <*.tmp> files.')

    additional = parser.add_argument_group('additional arguments - filtering')
    additional.add_argument(
        '-m',
        metavar='INT',
        type=int,
        default=25,
        help='Max. number of target sequences to report (--max-target-seqs/-k in diamond).')

    additional.add_argument(
        '-e',
        metavar='FLOAT',
        type=float,
        default=1e-5,
        help='Max. expected value to report alignments (--evalue/-e in diamond).')

    additional.add_argument(
        '-i',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Min. identity in percentage to report alignments. If "0" then set 90 - 2.5 * 100 * median sequence divergence.')

    additional.add_argument(
        '-s',
        metavar='FLOAT',
        type=float,
        default=90,
        help='Min. subject cover within a read cluster to report alignments.')

    additional.add_argument(
        '-n',
        metavar='INT',
        type=int,
        default=2147483647,
        help='Max. number of secondary alignments to report (-N in minimap2).')

    additional.add_argument(
        '-p',
        metavar='FLOAT',
        type=float,
        default=0.9,
        help='Min. secondary-to-primary score ratio to report secondary alignments (-p in minimap2).')

    additional.add_argument(
        '-z',
        metavar='FLOAT',
        type=float,
        default=1,
        help='Min. estimated genome copies of a species to report it ARG copies and abundances.')

    additional.add_argument(
        '-u',
        metavar='INT',
        type=int,
        default=0,
        help='Max. number of ARG-containing reads per chunk for overlapping. If "0" then use a single chunk.')

    additional = parser.add_argument_group('additional arguments - graph clustering')
    additional.add_argument(
        '-b',
        metavar='INT',
        type=int,
        default=1000,
        help='Terminal condition - max. iterations.')

    additional.add_argument(
        '-x',
        metavar='FLOAT',
        type=float,
        default=2,
        help='MCL parameter - inflation.')

    additional.add_argument(
        '-y',
        metavar='FLOAT',
        type=float,
        default=2,
        help='MCL parameter - expansion.')

    parser.add_argument('-v', '--version', action='version', version=__version__, help=argparse.SUPPRESS)
    parser.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)

    if len(argv) == 1:
        print(f"\
                       \n\
 ___ ________ ____     \n\
/ _ `/ __/ _ `/ _ \\   \n\
\\_,_/_/  \\_, /\\___/ \n\
        /___/       ver. {__version__}\n")

    opt = parser.parse_args(argv[1:])
    run(opt)


def run(opt):
    '''
    Sanity check of options.
    '''
    ## check for output folder
    if not os.path.isdir(opt.output):
        os.makedirs(opt.output, exist_ok=True)
    else:
        logger.warning(f'Folder <{opt.output}> exists. Files will be overwritten.')

    ## check for input files
    for file in opt.files:
        if not os.path.isfile(file):
            logger.critical(f'File <{file}> does not exist.')
            sys.exit(2)

    ## check for database
    if not os.path.isdir(opt.db):
        logger.critical(f'Database folder <{opt.db}> does not exist.')
        sys.exit(2)
    else:
        files = [os.path.basename(file) for file in glob.glob(os.path.join(opt.db, '*'))]
        if (
            'metadata.tsv' not in files or 'sarg.metadata.tsv' not in files
            or 'prot.dmnd' not in files or 'sarg.dmnd' not in files
            or len([file for file in files if 'nucl.' in file and '.mmi' in file]) != 16
            or len([file for file in files if 'sarg.' in file and '.mmi' in file]) <= 32
        ):
            logger.critical(f'Database <{opt.db}> is not complete or not indexed.')
            sys.exit(2)

    ## check for logical cores
    if opt.threads > os.cpu_count():
        logger.warning(
            f'Threads <{opt.threads}> exceeds available logical cores, will use <{os.cpu_count()}> instead.')
        opt.threads = os.cpu_count()

    ## run
    for index, file in enumerate(opt.files):
        if len(opt.files) > 1:
            logger.info(f'Processing file <{file}> ({index + 1}/{len(opt.files)}) ...')

        if not opt.skip_melon or any(not os.path.isfile(get_filename(file, opt.output, ext)) for ext in ['.tsv', '.json']):
            GenomeProfiler(file, opt.db, opt.output, opt.threads).run(skip_clean=opt.skip_clean)
        else:
            logger.info('Loading genome copy information ...')

        AntibioticResistanceGeneProfiler(file, opt.db, opt.output, opt.threads).run(
            plasmid=opt.plasmid, skip_clean=opt.skip_clean,
            max_target_seqs=opt.m, evalue=opt.e, identity=opt.i, subject_cover=opt.s,
            secondary_num=opt.n, secondary_ratio=opt.p,
            min_genome_copies=opt.z, chunk_size=opt.u,
            max_iterations=opt.b, inflation=opt.x, expansion=opt.y)

        if index == len(opt.files) - 1:
            logger.info('Done.')
        else:
            logger.info('Done.\n')


if __name__ == '__main__':
    cli(sys.argv)
