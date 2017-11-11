from inspect import stack
from os.path import dirname, isfile, join

from .process_bam import index_bam_using_samtools, sort_bam_using_samtools
from .support.support.path import remove_path
from .support.support.subprocess_ import run_command_and_monitor


def align_fastq_gzs_using_bwa(fasta_gz_file_path,
                              fastq_gz_file_paths,
                              n_jobs=1,
                              overwrite=False):
    """
    Align unpaired or paired .fastq.gz file using bwa.
    Arguments:
        fasta_gz_file_path (str):
        fastq_gz_file_paths (iterable): (<= 2) unpaired or paired sequences
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    if not all([
            isfile('{}.{}'.format(fasta_gz_file_path, suffix))
            for suffix in [
                'bwt',
                'pac',
                'ann',
                'amb',
                'sa',
            ]
    ]):
        print('Indexing ...')

        run_command_and_monitor(
            'bwa index {}'.format(fasta_gz_file_path), print_command=True)

    if 2 < len(fastq_gz_file_paths):
        raise ValueError(
            'fastq_gz_file_paths must contain unpaired or paired .fastq.gz file path.'
        )

    output_bam_file_path = join(
        dirname(fastq_gz_file_paths[0]), stack()[0][3] + '.bam')

    if isfile(output_bam_file_path) and not overwrite:
        raise FileExistsError('{} exists.'.format(output_bam_file_path))

    run_command_and_monitor(
        'bwa mem -t {} {} {} | samtools view -Sb --threads {} > {}'.format(
            n_jobs, fasta_gz_file_path, ' '.join(fastq_gz_file_paths), n_jobs,
            output_bam_file_path),
        print_command=True)

    return index_bam_using_samtools(
        sort_bam_using_samtools(
            output_bam_file_path, n_jobs=n_jobs, overwrite=overwrite),
        n_jobs=n_jobs,
        overwrite=overwrite)


def align_fastq_gzs_using_hisat2(fasta_file_path,
                                 fastq_gz_file_paths,
                                 sequence_type,
                                 n_jobs=1,
                                 overwrite=False):
    """
    Align unpaired or paired .fastq.gz files using hisat2.
    Arguments:
        fasta_file_path (str): reference .fasta.gz file path
        fastq_gz_file_paths (iterable): (<= 2) unpaired or paired end sequences
        sequence_type (str): 'DNA' | 'RNA'
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    if not all(
        [isfile('{}.{}.ht2'.format(fasta_file_path, i)) for i in range(1, 9)]):
        print('Indexing ...')

        run_command_and_monitor(
            'hisat2-build {0} {0}'.format(fasta_file_path), print_command=True)

    output_bam_file_path = join(
        dirname(fastq_gz_file_paths[0]), stack()[0][3] + '.bam')

    if isfile(output_bam_file_path) and not overwrite:
        raise FileExistsError('{} exists.'.format(output_bam_file_path))

    summary_file_path = output_bam_file_path + '.summary'

    if len(fastq_gz_file_paths) == 1:
        print('Using single-end ...')

        sample_argument = '-U {}'.format(*fastq_gz_file_paths)

    elif len(fastq_gz_file_paths) == 2:
        print('Using paired-end ...')

        sample_argument = '-1 {} -2 {}'.format(*fastq_gz_file_paths)

    else:
        raise ValueError(
            'fastq_gz_file_paths must contain unpaired or paired .fastq.gz file path.'
        )

    if sequence_type == 'DNA':
        additional_argument = '--no-spliced-alignment'

    elif sequence_type == 'RNA':
        additional_argument = '--dta --dta-cufflinks'

    run_command_and_monitor(
        'hisat2 {} -x {} --summary-file {} --threads {} {} | samtools view -Sb --threads {} > {}'.
        format(sample_argument, fasta_file_path, summary_file_path, n_jobs,
               additional_argument, n_jobs, output_bam_file_path),
        print_command=True)

    return index_bam_using_samtools(
        sort_bam_using_samtools(output_bam_file_path, n_jobs=n_jobs),
        n_jobs=n_jobs,
        overwrite=overwrite)


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

    if not isfile(fasta_gz_kallisto_index_file_path):
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

    if overwrite:
        remove_path(output_directory_path)

    run_command_and_monitor(
        'kallisto quant --index {} --output-dir {} --bootstrap-samples {} {}'.
        format(fasta_gz_kallisto_index_file_path, output_directory_path,
               n_bootstraps, sample_argument),
        print_command=True)

    return output_directory_path
