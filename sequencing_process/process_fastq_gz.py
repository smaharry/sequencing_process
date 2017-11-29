from inspect import stack
from os.path import dirname, exists, join
from sys import platform

from . import RESOURCE_DIRECTORY_PATH
from .support.support.subprocess_ import run_command


def check_fastq_gzs_using_fastqc(fastq_gz_file_paths,
                                 n_jobs=1,
                                 overwrite=False):
    """
    Check .fastq.gz files using fastqc.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) 1 (unpaired) or 2 (paired)
            .fastq.gz file path
        n_jobs (int):
        overwrite (bool):
    Returns:
        None
    """

    for fp in fastq_gz_file_paths:
        fp += '_fastqc.html'
        if not overwrite and exists(fp):
            raise FileExistsError(fp)

    check_fastq_gzs(fastq_gz_file_paths)

    run_command(
        'fastqc --threads {} {}'.format(n_jobs, ' '.join(fastq_gz_file_paths)),
        print_command=True)


def trim_fastq_gzs_using_skewer(fastq_gz_file_paths,
                                forward_bad_sequence_fasta_file_path=join(
                                    RESOURCE_DIRECTORY_PATH,
                                    'general_bad_sequence.fasta'),
                                reverse_bad_sequence_fasta_file_path=join(
                                    RESOURCE_DIRECTORY_PATH,
                                    'general_bad_sequence.fasta'),
                                snv_error_rate=0,
                                indel_error_rate=0,
                                overlap_length=12,
                                end_quality=30,
                                min_length_after_trimming=30,
                                remove_n=True,
                                n_jobs=1,
                                output_prefix=None,
                                overwrite=False):
    """
    Trim .fastq.gz files using skewer.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) 1 (unpaired) or 2 (paired)
            .fastq.gz file path
        forward_bad_sequence_fasta_file_path (str):
        reverse_bad_sequence_fasta_file_path (str):
        snv_error_rate (float):
        indel_error_rate (float):
        overlap_length (int):
        end_quality (number):
        min_length_after_trimming (int):
        remove_n (bool):
        n_jobs (int):
        output_prefix (str):
        overwrite (bool):
    Returns:
        list:
    """

    check_fastq_gzs(fastq_gz_file_paths)

    if not output_prefix:
        output_prefix = join(dirname(fastq_gz_file_paths[0]), stack()[0][3])

    output_fastq_file_paths = [
        '{}-trimmed-pair{}.fastq'.format(output_prefix, i) for i in [1, 2]
    ]
    for fp in output_fastq_file_paths:
        if not overwrite and exists(fp + '.gz'):
            raise FileExistsError(fp + '.gz')

    command = 'skewer -x {} -r {} -d {} --end-quality {} --min {} {} --output {} --masked-output --excluded-output --threads {}'.format(
        forward_bad_sequence_fasta_file_path, snv_error_rate, indel_error_rate,
        end_quality, min_length_after_trimming, ['', '-n'][remove_n],
        output_prefix, n_jobs)

    additional_arguments = []
    if len(fastq_gz_file_paths) == 1:
        additional_arguments.append('-k {}'.format(overlap_length))
    else:
        additional_arguments.append(
            '-y {}'.format(reverse_bad_sequence_fasta_file_path))
    additional_arguments.extend(fastq_gz_file_paths)

    run_command(
        '{} {}'.format(command, ' '.join(additional_arguments)),
        print_command=True)

    log_file_path = output_prefix + '-trimmed.log'
    print('{}:'.format(log_file_path))
    with open(log_file_path) as f:
        print(f.read())

    for fp in output_fastq_file_paths:
        run_command('gzip --force {}'.format(fp), print_command=True)

    return [fp + '.gz' for fp in output_fastq_file_paths]


def align_fastq_gzs_using_bwa_mem(fastq_gz_file_paths,
                                  fasta_gz_file_path,
                                  n_jobs=1,
                                  output_bam_file_path=None,
                                  overwrite=False):
    """
    Align .fastq.gz files using bwa mem.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) 1 (unpaired) or 2 (paired)
            .fastq.gz file path
        fasta_gz_file_path (str):
        n_jobs (int):
        output_bam_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    check_fastq_gzs(fastq_gz_file_paths)

    if not all([
            exists(fasta_gz_file_path + e)
            for e in ['.bwt', '.pac', '.ann', '.amb', '.sa']
    ]):
        run_command(
            'bwa index {}'.format(fasta_gz_file_path), print_command=True)

    if not exists(fasta_gz_file_path + '.alt'):
        raise FileNotFoundError('ALT-aware BWA-MEM alignment needs {}.'.format(
            fasta_gz_file_path + '.alt'))

    if not output_bam_file_path:
        output_bam_file_path = join(
            dirname(fastq_gz_file_paths[0]), stack()[0][3] + '.bam')

    if not overwrite and exists(output_bam_file_path):
        raise FileExistsError(output_bam_file_path)

    run_command(
        'bwa mem -t {} {} {} | {} {} {}.alt | samtools view -Sb --threads {} > {}'.
        format(n_jobs, fasta_gz_file_path, ' '.join(fastq_gz_file_paths),
               join(RESOURCE_DIRECTORY_PATH, 'k8-0.2.3',
                    'k8-{}'.format(platform)),
               join(RESOURCE_DIRECTORY_PATH, 'bwa-postalt.js'),
               fasta_gz_file_path, n_jobs, output_bam_file_path),
        print_command=True)

    return output_bam_file_path


def align_fastq_gzs_using_hisat2(fastq_gz_file_paths,
                                 fasta_file_path,
                                 sequence_type,
                                 n_jobs=1,
                                 output_bam_file_path=None,
                                 overwrite=False):
    """
    Align .fastq.gz files using hisat2.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) 1 (unpaired) or 2 (paired)
            .fastq.gz file path
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

    check_fastq_gzs(fastq_gz_file_paths)

    if not all([
            exists(fasta_file_path + '.{}.ht2'.format(i))
            for i in [1, 2, 3, 4, 5, 6, 7, 8]
    ]):
        run_command(
            'hisat2-build {0} {0}'.format(fasta_file_path), print_command=True)

    if not output_bam_file_path:
        output_bam_file_path = join(
            dirname(fastq_gz_file_paths[0]), stack()[0][3] + '.bam')

    if not overwrite and exists(output_bam_file_path):
        raise FileExistsError(output_bam_file_path)

    additional_arguments = []
    if len(fastq_gz_file_paths) == 1:
        additional_arguments.append('-U {}'.format(*fastq_gz_file_paths))
    else:
        additional_arguments.append('-1 {} -2 {}'.format(*fastq_gz_file_paths))
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


def count_transcripts_using_kallisto_quant(
        fastq_gz_file_paths,
        fasta_gz_file_path,
        output_directory_path,
        n_bootstraps=100,
        fragment_lendth=180,
        fragment_lendth_standard_deviation=20,
        n_jobs=1,
        overwrite=False):
    """
    Count transcripts using kallisto quant.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) 1 (unpaired) or 2 (paired)
            .fastq.gz file path
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

    check_fastq_gzs(fastq_gz_file_paths)

    fasta_gz_kallisto_index_file_path = '{}.kallisto.index'.format(
        fasta_gz_file_path)
    if not exists(fasta_gz_kallisto_index_file_path):
        run_command(
            'kallisto index --index {} {}'.format(
                fasta_gz_kallisto_index_file_path, fasta_gz_file_path),
            print_command=True)

    if len(fastq_gz_file_paths) == 1:
        sample_argument = '--single --fragment-length {} --sd {} {}'.format(
            fragment_lendth, fragment_lendth_standard_deviation,
            *fastq_gz_file_paths)
    else:
        sample_argument = '{} {}'.format(*fastq_gz_file_paths)

    if not overwrite and exists(output_directory_path):
        raise FileExistsError(output_directory_path)

    run_command(
        'kallisto quant --index {} --output-dir {} --bootstrap-samples {} {}'.
        format(fasta_gz_kallisto_index_file_path, output_directory_path,
               n_bootstraps, sample_argument),
        print_command=True)

    return output_directory_path


def check_fastq_gzs(fastq_gz_file_paths):
    """
    Check .fastq.gz files.
    Arguments:
        fastq_gz_file_paths (iterable): (<= 2) 1 (unpaired) or 2 (paired)
            .fastq.gz file path
    Returns:
        None
    """

    if len(fastq_gz_file_paths) not in [1, 2]:
        raise ValueError(
            'fastq_gz_file_paths must contain either 1 (unpaired) or 2 (paired) .fastq.gz file path.'
        )
