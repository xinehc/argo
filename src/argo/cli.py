import sys
import os
import glob

from argparse import ArgumentParser, SUPPRESS
from . import __version__
from melon.utils import logger
from melon import GenomeProfiler
from . import AntibioticResistanceGeneProfiler


def cli(argv=sys.argv):
    '''
    Entry point for command line interface.
    '''
    parser = ArgumentParser(description='ARGO: Overlapping-based ARG profiling', add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    additional = parser.add_argument_group('additional arguments for profiling genomes')
    additional_arg = parser.add_argument_group('additional arguments for profiling antibiotic resistance genes')

    parser.add_argument(
        'FILE',
        nargs='+',
        help='Input fasta <*.fa|*.fasta> or fastq <*.fq|*.fastq> file, gzip optional <*.gz>.')

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

    optional.add_argument(
        '-t',
        '--threads',
        metavar='INT',
        type=int,
        default=os.cpu_count(),
        help=f'Number of threads. [{os.cpu_count()}]')

    optional.add_argument(
        '-k',
        '--db-kraken',
        metavar='DIR',
        help='Unzipped kraken2 database for pre-filtering of non-prokaryotic reads. Skip if not given.')

    optional.add_argument(
        '--skip-clean',
        action='store_true',
        help='Skip cleaning, keep all temporary <*.tmp> files.')

    additional.add_argument(
        '-m',
        metavar='INT',
        type=int,
        default=25,
        help='Max. number of target sequences to report (--max-target-seqs/-k in diamond). [25]')

    additional.add_argument(
        '-e',
        metavar='FLOAT',
        type=float,
        default=1e-15,
        help='Max. expected value to report alignments (--evalue/-e in diamond). [1e-15]')

    additional.add_argument(
        '-i',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Min. identity in percentage to report alignments (--id in diamond). [0]')

    additional.add_argument(
        '-s',
        metavar='FLOAT',
        type=float,
        default=75,
        help='Min. subject cover to report alignments (--subject-cover in diamond). [75]')

    additional.add_argument(
        '-n',
        metavar='INT',
        type=int,
        default=2147483647,
        help='Max. number of secondary alignments to report (-N in minimap2). [2147483647]')

    additional.add_argument(
        '-p',
        metavar='FLOAT',
        type=float,
        default=0.9,
        help='Min. secondary-to-primary score ratio to report secondary alignments (-p in minimap2). [0.9]')

    additional.add_argument(
        '-a',
        metavar='INT',
        type=int,
        default=1000,
        help='Terminal condition for EM - max. iterations. [1000]')

    additional.add_argument(
        '-c',
        metavar='FLOAT',
        type=float,
        default=1e-10,
        help='Terminal condition for EM - epsilon (precision). [1e-10]')

    additional_arg.add_argument(
        '-M',
        metavar='INT',
        type=int,
        default=25,
        help='Max. number of target sequences to report (--max-target-seqs/-k in diamond). [25]')

    additional_arg.add_argument(
        '-E',
        metavar='FLOAT',
        type=float,
        default=1e-5,
        help='Max. expected value to report alignments (--evalue/-e in diamond). [1e-5]')

    additional_arg.add_argument(
        '-I',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Min. identity in percentage to report alignments. If 0 then 90 - 2.5 * median sequence divergence in percentage. [0]')

    additional_arg.add_argument(
        '-S',
        metavar='FLOAT',
        type=float,
        default=90,
        help='Min. subject cover of all HPSs within a read cluster to report alignments. [90]')

    additional_arg.add_argument(
        '-N',
        metavar='INT',
        type=int,
        default=2147483647,
        help='Max. number of secondary alignments to report (-N in minimap2). [2147483647]')

    additional_arg.add_argument(
        '-P',
        metavar='FLOAT',
        type=float,
        default=0.9,
        help='Min. secondary-to-primary score ratio to report secondary alignments (-p in minimap2). [0.9]')

    additional_arg.add_argument(
        '-z',
        metavar='FLOAT',
        type=float,
        default=1,
        help='Min. estimated genome copies of a species to report it ARG copies and abundance. [1]')

    additional_arg.add_argument(
        '-u',
        metavar='INT',
        type=int,
        default=0,
        help='Max. number of ARG-containing reads per chunk for overlapping. If 0 then use a single chunk. [0]')

    additional_arg.add_argument(
        '-b',
        metavar='INT',
        type=int,
        default=1000,
        help='Graph clustering parameter for MCL - max. iterations. [1000]')

    additional_arg.add_argument(
        '-x',
        metavar='FLOAT',
        type=float,
        default=2,
        help='Graph clustering parameter for MCL - inflation. [2]')

    additional_arg.add_argument(
        '-y',
        metavar='FLOAT',
        type=float,
        default=2,
        help='Graph clustering parameter for MCL - expansion. [2]')

    parser.add_argument('-v', '--version', action='version', version=__version__, help=SUPPRESS)
    parser.add_argument('-h', '--help', action='help', help=SUPPRESS)

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
    for file in opt.FILE:
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
            or 'prot.dmnd' not in files  or 'sarg.dmnd' not in files 
            or len([file for file in files if 'nucl' in file and '.mmi' in file]) != 16
            or len([file for file in files if 'sarg' in file and '.mmi' in file]) != 43
        ):
            logger.critical(f'Database <{opt.db}> is not complete or not indexed.')
            sys.exit(2)

    ## check for kraken2 database
    if opt.db_kraken is not None:
        if not os.path.isdir(opt.db_kraken):
            logger.critical(f'Kraken2 database folder <{opt.db_kraken}> does not exist.')
            sys.exit(2)
        else:
            files = [os.path.basename(file) for file in glob.glob(os.path.join(opt.db_kraken, '*'))]
            if 'ktaxonomy.tsv' not in files or len([file for file in files if 'database' in file]) != 7:
                logger.critical(f'Kraken2 database <{opt.db_kraken}> is not complete.')
                sys.exit(2)

    ## check for logical cores
    if opt.threads > os.cpu_count():
        logger.warning(
            f'Threads <{opt.threads}> exceeds available logical cores, will use <{os.cpu_count()}> instead.')
        opt.threads = os.cpu_count()
    os.environ['OMP_NUM_THREADS'] = str(opt.threads)

    ## run
    for index, file in enumerate(opt.FILE):
        if len(opt.FILE) > 1:
            logger.info(f'Processing file <{file}> ({index + 1}/{len(opt.FILE)}) ...')

        GenomeProfiler(file, opt.db, opt.output, opt.threads).run(
            db_kraken=opt.db_kraken, skip_clean=opt.skip_clean,
            max_target_seqs=opt.m, evalue=opt.e, identity=opt.i, subject_cover=opt.s,
            secondary_num=opt.n, secondary_ratio=opt.p,
            max_iterations=opt.a, epsilon=opt.c)

        AntibioticResistanceGeneProfiler(file, opt.db, opt.output, opt.threads).run(
            skip_clean=opt.skip_clean,
            max_target_seqs=opt.M, evalue=opt.E, identity=opt.I, subject_cover=opt.S,
            secondary_num=opt.N, secondary_ratio=opt.P,
            min_genome_copies=opt.z, chunk_size=opt.u,
            max_iterations=opt.b, inflation=opt.x, expansion=opt.y)

        if index == len(opt.FILE) - 1:
            logger.info('Done.')
        else:
            logger.info('Done.\n')


if __name__ == '__main__':
    cli(sys.argv)
