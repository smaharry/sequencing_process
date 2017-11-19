from inspect import stack
from os.path import dirname, exists, join
from sys import platform

from .support.support.path import clean_path
from .support.support.subprocess_ import run_command, run_command_and_monitor


def validate_fastq_gz_using_fqtools(fastq_gz_file_path):
    """
    Validate .fastq.gz file using fqtools.
    Arguments:
        fastq_gz_file_path (str):
    Returns:
        bool:
    """

    return run_command(
        'fqtools validate {}'.format(fastq_gz_file_path),
        print_command=True).stdout.strip()


def align_fastq_gzs_using_bwa(fasta_gz_file_path,
                              fastq_gz_file_paths,
                              n_jobs=1,
                              output_bam_file_path=None,
                              overwrite=False):
    """
    Align unpaired or paired .fastq.gz file using bwa.
    Arguments:
        fasta_gz_file_path (str):
        fastq_gz_file_paths (iterable): (<= 2) unpaired or paired sequences
        n_jobs (int):
        output_bam_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not all([
            exists('{}.{}'.format(fasta_gz_file_path, suffix))
            for suffix in ['bwt', 'pac', 'ann', 'amb', 'sa']
    ]):
        print('Indexing ...')
        run_command_and_monitor(
            'bwa index {}'.format(fasta_gz_file_path), print_command=True)

    if not exists('{}.alt'.format(fasta_gz_file_path)):
        raise FileExistsError('ALT-aware BWA-MEM alignment needs {}.'.format(
            output_bam_file_path))

    if 2 < len(fastq_gz_file_paths):
        raise ValueError(
            'fastq_gz_file_paths must contain unpaired or paired .fastq.gz file path.'
        )

    if not output_bam_file_path:
        output_bam_file_path = join(
            dirname(fastq_gz_file_paths[0]), stack()[0][3] + '.bam')

    if not overwrite and exists(output_bam_file_path):
        raise FileExistsError(output_bam_file_path)

    directory_path = dirname(clean_path(__file__))
    k8_path = join(directory_path, 'k8-0.2.3', 'k8-{}'.format(platform))
    js_path = join(directory_path, 'bwa-postalt.js')

    run_command_and_monitor(
        'bwa mem -t {} {} {} | {} {} {}.alt | samtools view -Sb --threads {} > {}'.
        format(n_jobs, fasta_gz_file_path, ' '.join(fastq_gz_file_paths),
               k8_path, js_path, fasta_gz_file_path, n_jobs,
               output_bam_file_path),
        print_command=True)

    return output_bam_file_path


def align_fastq_gzs_using_hisat2(fasta_file_path,
                                 fastq_gz_file_paths,
                                 sequence_type,
                                 n_jobs=1,
                                 output_bam_file_path=None,
                                 overwrite=False):
    """
    Align unpaired or paired .fastq.gz files using hisat2.
    Arguments:
        fasta_file_path (str): reference .fasta.gz file path
        fastq_gz_file_paths (iterable): (<= 2) unpaired or paired end sequences
        sequence_type (str): 'DNA' | 'RNA'
        n_jobs (int):
        output_bam_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not all(
        [exists('{}.{}.ht2'.format(fasta_file_path, i)) for i in range(1, 9)]):
        print('Indexing ...')
        run_command_and_monitor(
            'hisat2-build {0} {0}'.format(fasta_file_path), print_command=True)

    if len(fastq_gz_file_paths) == 2:
        raise ValueError(
            'fastq_gz_file_paths must contain paired .fastq.gz file paths.')

    if not output_bam_file_path:
        output_bam_file_path = join(
            dirname(fastq_gz_file_paths[0]), stack()[0][3] + '.bam')

    if not overwrite and exists(output_bam_file_path):
        raise FileExistsError(output_bam_file_path)

    additional_arguments = []

    if len(fastq_gz_file_paths) == 1:
        print('Using single-end ...')
        additional_arguments.append('-U {}'.format(*fastq_gz_file_paths))

    elif len(fastq_gz_file_paths) == 2:
        print('Using paired-end ...')
        additional_arguments.append('-1 {} -2 {}'.format(*fastq_gz_file_paths))

    else:
        raise ValueError(
            'fastq_gz_file_paths must contain unpaired or paired .fastq.gz file path.'
        )

    if sequence_type == 'DNA':
        additional_arguments.append('--no-spliced-alignment')
    elif sequence_type == 'RNA':
        additional_arguments.append('--dta --dta-cufflinks')
    else:
        raise ValueError('Unknown sequence_type {}.'.format(sequence_type))

    run_command_and_monitor(
        'hisat2 -x {} --summary-file {}.summary --threads {} {} | samtools view -Sb --threads {} > {}'.
        format(fasta_file_path, output_bam_file_path, n_jobs,
               ' '.join(additional_arguments), n_jobs, output_bam_file_path),
        print_command=True)

    return output_bam_file_path


def count_transcripts_using_kallisto(fasta_gz_file_path,
                                     fastq_gz_file_paths,
                                     output_directory_path,
                                     n_bootstraps=100,
                                     fragment_lendth=180,
                                     fragment_lendth_standard_deviation=20,
                                     n_jobs=1,
                                     overwrite=False):
    """
    Count transcripts using kallisto.
    Arguments:
        fasta_gz_file_path (str): cDNA sequences
        fastq_gz_file_paths (iterable): (<= 2)
        output_directory_path (str):
        n_bootstraps (int):
        fragment_lendth (number): estimated fragment length
        fragment_lendth_standard_deviation (number): estimated fragment length
            standard deviation
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    fasta_gz_kallisto_index_file_path = '{}.kallisto.index'.format(
        fasta_gz_file_path)

    if not exists(fasta_gz_kallisto_index_file_path):
        print('Indexing ...')
        run_command_and_monitor(
            'kallisto index --index {} {}'.format(
                fasta_gz_kallisto_index_file_path, fasta_gz_file_path),
            print_command=True)

    if len(fastq_gz_file_paths) == 1:
        print('Using single-end ...')
        sample_argument = '--single --fragment-length {} --sd {} {}'.format(
            fragment_lendth, fragment_lendth_standard_deviation,
            *fastq_gz_file_paths)

    elif len(fastq_gz_file_paths) == 2:
        print('Using paired-end ...')
        sample_argument = '{} {}'.format(*fastq_gz_file_paths)

    else:
        raise ValueError(
            'fastq_gz_file_paths must contain unpaired or paired .fastq.gz file path.'
        )

    if not overwrite and exists(output_directory_path):
        raise FileExistsError(output_directory_path)

    run_command_and_monitor(
        'kallisto quant --index {} --output-dir {} --bootstrap-samples {} {}'.
        format(fasta_gz_kallisto_index_file_path, output_directory_path,
               n_bootstraps, sample_argument),
        print_command=True)

    return output_directory_path
