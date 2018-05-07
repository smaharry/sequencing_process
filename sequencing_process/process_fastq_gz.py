from inspect import stack
from os.path import dirname, isdir, isfile, join
from sys import platform

from . import RESOURCE_DIRECTORY_PATH
from .print_and_run_command import print_and_run_command
from .support.support.multiprocess import multiprocess


def check_fastq_gzs(fastq_gz_file_paths):

    if len(fastq_gz_file_paths) not in (1, 2):

        raise ValueError(
            'fastq_gz_file_paths should contain either 1 (unpaired) or 2 (paired) .fastq.gz file paths.'
        )


def check_fastq_gzs_using_fastqc(fastq_gz_file_paths, n_job=1,
                                 overwrite=False):

    for fastq_gz_file_path in fastq_gz_file_paths:

        html_file_path = fastq_gz_file_path + '_fastqc.html'

        if not overwrite and isfile(html_file_path):

            raise FileExistsError(html_file_path)

    print_and_run_command('fastqc --threads {} {}'.format(
        n_job, ' '.join(fastq_gz_file_paths)))


def trim_fastq_gzs_using_skewer(fastq_gz_file_paths,
                                forward_bad_sequences_fasta_file_path=join(
                                    RESOURCE_DIRECTORY_PATH,
                                    'general_bad_sequences.fasta'),
                                reverse_bad_sequences_fasta_file_path=join(
                                    RESOURCE_DIRECTORY_PATH,
                                    'general_bad_sequences.fasta'),
                                snv_error_rate=0,
                                indel_error_rate=0,
                                overlap_length=12,
                                end_quality=30,
                                min_length_after_trimming=30,
                                remove_n=True,
                                n_job=1,
                                output_directory_path=None,
                                overwrite=False):

    check_fastq_gzs(fastq_gz_file_paths)

    if output_directory_path is None:

        output_directory_path = join(
            dirname(fastq_gz_file_paths[0]),
            stack()[0][3])

    if not output_directory_path.endswith('/'):

        output_directory_path += '/'

    if not overwrite and isdir(output_directory_path):

        raise FileExistsError(output_directory_path)

    additional_arguments = []

    if len(fastq_gz_file_paths) == 1:

        additional_arguments.append('-m tail')

        additional_arguments.append('-k {}'.format(overlap_length))

    elif len(fastq_gz_file_paths) == 2:

        additional_arguments.append('-m pe')

        additional_arguments.append(
            '-y {}'.format(reverse_bad_sequences_fasta_file_path))

    print_and_run_command(
        'skewer -x {} -r {} -d {} --end-quality {} --min {} {} --output {} --masked-output --excluded-output --threads {} {}'.
        format(forward_bad_sequences_fasta_file_path, snv_error_rate,
               indel_error_rate, end_quality, min_length_after_trimming,
               ('', '-n')[remove_n], output_directory_path, n_job,
               ' '.join(additional_arguments + list(fastq_gz_file_paths))))

    log_file_path = join(output_directory_path, 'trimmed.log')

    print('{}:'.format(log_file_path))

    with open(log_file_path) as log_file:

        print(log_file.read())

    output_fastq_file_paths = [
        join(output_directory_path, 'trimmed-pair{}.fastq'.format(i))
        for i in (1, 2)
    ]

    return multiprocess(
        _gzip_compress,
        ((outptu_fastq_file_path, )
         for outptu_fastq_file_path in output_fastq_file_paths),
        n_job=n_job)


def _gzip_compress(file_path):

    print_and_run_command('gzip --force {}'.format(file_path))

    return file_path + '.gz'


def align_fastq_gzs_using_bwa_mem(fastq_gz_file_paths,
                                  fasta_gz_file_path,
                                  n_job=1,
                                  output_bam_file_path=None,
                                  overwrite=False):

    check_fastq_gzs(fastq_gz_file_paths)

    if not all(
            isfile(fasta_gz_file_path + extension)
            for extension in ('.bwt', '.pac', '.ann', '.amb', '.sa')):

        print_and_run_command('bwa index {}'.format(fasta_gz_file_path))

    if not isfile(fasta_gz_file_path + '.alt'):

        raise FileNotFoundError('ALT-aware BWA-MEM alignment needs {}.'.format(
            fasta_gz_file_path + '.alt'))

    if output_bam_file_path is None:

        output_bam_file_path = join(
            dirname(fastq_gz_file_paths[0]),
            stack()[0][3] + '.bam')

    if not overwrite and isfile(output_bam_file_path):

        raise FileExistsError(output_bam_file_path)

    print_and_run_command(
        'bwa mem -t {} -v 3 {} {} | {} {} {}.alt | samtools view -Sb --threads {} > {}'.
        format(n_job, fasta_gz_file_path, ' '.join(fastq_gz_file_paths),
               join(RESOURCE_DIRECTORY_PATH, 'k8-0.2.3',
                    'k8-{}'.format(platform)),
               join(RESOURCE_DIRECTORY_PATH, 'bwa-postalt.js'),
               fasta_gz_file_path, n_job, output_bam_file_path))

    return output_bam_file_path


def align_fastq_gzs_using_hisat2(fastq_gz_file_paths,
                                 fasta_file_path,
                                 sequence_type,
                                 n_job=1,
                                 output_bam_file_path=None,
                                 overwrite=False):

    check_fastq_gzs(fastq_gz_file_paths)

    if not all(
            isfile(fasta_file_path + '.{}.ht2'.format(i))
            for i in (1, 2, 3, 4, 5, 6, 7, 8)):

        print_and_run_command('hisat2-build {0} {0}'.format(fasta_file_path))

    if output_bam_file_path is None:

        output_bam_file_path = join(
            dirname(fastq_gz_file_paths[0]),
            stack()[0][3] + '.bam')

    if not overwrite and isfile(output_bam_file_path):

        raise FileExistsError(output_bam_file_path)

    additional_arguments = []

    if len(fastq_gz_file_paths) == 1:

        additional_arguments.append('-U {}'.format(fastq_gz_file_paths[0]))

    elif len(fastq_gz_file_paths) == 2:

        additional_arguments.append('-1 {} -2 {}'.format(*fastq_gz_file_paths))

    if sequence_type not in ('DNA', 'RNA'):

        raise ValueError('Unknown sequence_type: {}.'.format(sequence_type))

    elif sequence_type == 'DNA':

        additional_arguments.append('--no-spliced-alignment')

    elif sequence_type == 'RNA':

        additional_arguments.append('--dta --dta-cufflinks')

    print_and_run_command(
        'hisat2 -x {} --summary-file {} --threads {} {} | samtools view -Sb --threads {} > {}'.
        format(fasta_file_path, output_bam_file_path + '.summary', n_job,
               ' '.join(additional_arguments), n_job, output_bam_file_path))

    return output_bam_file_path


def count_transcripts_using_kallisto_quant(
        fastq_gz_file_paths,
        fasta_gz_file_path,
        output_directory_path,
        n_bootstrap=100,
        fragment_length=180,
        fragment_length_standard_deviation=20,
        n_job=1,
        overwrite=False):

    check_fastq_gzs(fastq_gz_file_paths)

    fasta_gz_kallisto_index_file_path = '{}.kallisto.index'.format(
        fasta_gz_file_path)

    if not isfile(fasta_gz_kallisto_index_file_path):

        print_and_run_command('kallisto index --index {} {}'.format(
            fasta_gz_kallisto_index_file_path, fasta_gz_file_path))

    if not overwrite and isdir(output_directory_path):

        raise FileExistsError(output_directory_path)

    if len(fastq_gz_file_paths) == 1:

        sample_argument = '--single --fragment-length {} --sd {} {}'.format(
            fragment_length, fragment_length_standard_deviation,
            fastq_gz_file_paths[0])

    elif len(fastq_gz_file_paths) == 2:

        sample_argument = '{} {}'.format(*fastq_gz_file_paths)

    print_and_run_command(
        'kallisto quant --index {} --output-dir {} --bootstrap-samples {} --threads {} {}'.
        format(fasta_gz_kallisto_index_file_path, output_directory_path,
               n_bootstrap, n_job, sample_argument))

    return output_directory_path
