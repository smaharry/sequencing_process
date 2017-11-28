from inspect import stack
from os.path import dirname, exists, join
from sys import platform

from .support.support.path import clean_path
from .support.support.subprocess_ import run_command


def trim_fastq_gzs_using_skewer(fastq_gz_file_paths,
                                bad_sequence_fasta_file_path=join(
                                    dirname(dirname(__file__)), 'resource',
                                    'general_bad_sequence.fasta'),
                                snv_error_rate=0,
                                indel_error_rate=0,
                                overlap_length=13,
                                end_quality=30,
                                min_length_after_trimming=30,
                                remove_n=True,
                                n_jobs=1,
                                output_fastq_gz_file_path_prefix=None):
    """
    Trim paired .fastq.gz files using skewer.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) unpaired or paired sequences
        bad_sequence_fasta_file_path (str):
        snv_error_rate (float):
        indel_error_rate (float):
        overlap_length (int):
        end_quality (number):
        min_length_after_trimming (int):
        remove_n (bool):
        n_jobs (int):
        output_fastq_gz_file_path_prefix (str):
    Returns:
        list:
    """

    check_fastq_gz_file_paths(fastq_gz_file_paths)

    common_command = 'skewer -x {} -y {} -r {} -d {} --end-quality {} --min {} {} --output {} --masked-output --excluded-output --threads {}'.
        format(bad_sequence_fasta_file_path, bad_sequence_fasta_file_path,
               snv_error_rate, indel_error_rate, end_quality,
               min_length_after_trimming, ['', '-n'][remove_n],
               output_fastq_gz_file_path_prefix, n_jobs)

    additional_arguments = []
    if len(fastq_gz_file_paths) == 1:
        additional_arguments.append('-k {}'.format(overlap_length))
    additional_arguments.extend(fastq_gz_file_paths)

    run_command('{} {}'.format(
        common_command, ' '.join(additional_arguments)),
        print_command=True)

    output_fastq_file_paths = [
        '{}-trimmed-pair{}.fastq'.format(output_fastq_gz_file_path_prefix, i)
        for i in [1, 2]
    ]

    for fp in output_fastq_file_paths:
        run_command('gzip {}'.format(fp), print_command=True)

    return ['{}.gz'.format(fp) for fp in output_fastq_file_paths]


def align_fastq_gzs_using_bwa(fastq_gz_file_paths,
                              fasta_gz_file_path,
                              n_jobs=1,
                              output_bam_file_path=None,
                              overwrite=False):
    """
    Align unpaired or paired .fastq.gz file using bwa.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) unpaired or paired sequences
        fasta_gz_file_path (str):
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
        run_command(
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

    run_command(
        'bwa mem -t {} {} {} | {} {} {}.alt | samtools view -Sb --threads {} > {}'.
        format(n_jobs, fasta_gz_file_path, ' '.join(fastq_gz_file_paths),
               k8_path, js_path, fasta_gz_file_path, n_jobs,
               output_bam_file_path),
        print_command=True)

    return output_bam_file_path


def align_fastq_gzs_using_hisat2(fastq_gz_file_paths,
                                 fasta_file_path,
                                 sequence_type,
                                 n_jobs=1,
                                 output_bam_file_path=None,
                                 overwrite=False):
    """
    Align unpaired or paired .fastq.gz files using hisat2.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) unpaired or paired sequences
        fasta_gz_file_path (str):
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
        run_command(
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

    run_command(
        'hisat2 -x {} --summary-file {}.summary --threads {} {} | samtools view -Sb --threads {} > {}'.
        format(fasta_file_path, output_bam_file_path, n_jobs,
               ' '.join(additional_arguments), n_jobs, output_bam_file_path),
        print_command=True)

    return output_bam_file_path


def count_transcripts_using_kallisto(fastq_gz_file_paths,
                                     fasta_gz_file_path,
                                     output_directory_path,
                                     n_bootstraps=100,
                                     fragment_lendth=180,
                                     fragment_lendth_standard_deviation=20,
                                     n_jobs=1,
                                     overwrite=False):
    """
    Count transcripts using kallisto.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) unpaired or paired sequences
        fasta_gz_file_path (str): cDNA sequences
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
        run_command(
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

    run_command(
        'kallisto quant --index {} --output-dir {} --bootstrap-samples {} {}'.
        format(fasta_gz_kallisto_index_file_path, output_directory_path,
               n_bootstraps, sample_argument),
        print_command=True)

    return output_directory_path


def check_fastq_gz_file_paths(fastq_gz_file_paths):
    """
    Check .fastq_gz_file_paths.
    Arguments:
        fastq_gz_file_paths (iterable):
    Returns:
        None
    """

    if len(fastq_gz_file_paths) == 1:
        print('Using single-end ...')

    elif len(fastq_gz_file_paths) == 2:
        print('Using paired-end ...')

    else:
        raise ValueError(
            'fastq_gz_file_paths must contain either 1 (unpaired) or 2 (paired) .fastq.gz file path.'
        )
