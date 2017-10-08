from os.path import exists

from .process_bam import index_bam_using_samtools, sort_bam_using_samtools
from .support.support.subprocess_ import run_command


def align_fastq_gzs_using_hisat2(hisat2_index_file_paths_prefix,
                                 fastq_gz_file_paths,
                                 sequence_type,
                                 n_jobs=1):
    """
    Align .fastq.gz files using hisat2.
    Arguments:
        hisat2_index_file_paths_prefix (str): reference .fasta file path
        fastq_gz_file_paths (iterable): (< 2)
        sequence_type (str): 'DNA' | 'RNA'
        n_jobs (int):
    Returns:
        str:
    """

    if not all([
            exists('{}.{}.ht2'.format(hisat2_index_file_paths_prefix, i))
            for i in range(1, 9)
    ]):
        print('Could not find HISAT2 index; building it ...')
        run_command('hisat2-build {} {}'.format(
            * [hisat2_index_file_paths_prefix] * 2))

    if len(fastq_gz_file_paths) == 1:
        print('Using single-end ...')
        sample_command = '-U {}'.format(*fastq_gz_file_paths)

    elif len(fastq_gz_file_paths) == 2:
        print('Using paired-end ...')
        sample_command = '-1 {} -2 {}'.format(*fastq_gz_file_paths)

    else:
        raise ValueError(
            'fastq_gz_file_paths must contain unpaired or paired .fastq.gz file paths.'
        )

    if sequence_type == 'DNA':
        additional_arguments = '--no-spliced-alignment'
    elif sequence_type == 'RNA':
        additional_arguments = '--dta --dta-cufflinks'

    output_bam_file_path = fastq_gz_file_paths[0] + '.align_fastq_gzs_using_hisat2.bam'

    run_command('hisat2 {} -x {} -p {} {} | samtools view -Sb -@ {} > {}'.
                format(sample_command, hisat2_index_file_paths_prefix, n_jobs,
                       additional_arguments, n_jobs, output_bam_file_path))

    return index_bam_using_samtools(
        sort_bam_using_samtools(output_bam_file_path, n_jobs=n_jobs),
        n_jobs=n_jobs)


def align_fastq_gzs_using_bwa(fasta_file_path, fastq_gz_file_paths, n_jobs=1):
    """
    Align unpaired or paired .fastq.gz files using bwa.
    Arguments:
        fasta_file_path (str):
        fastq_gz_file_paths (iterable): (< 2)
        n_jobs (int):
    Returns:
        str:
    """

    if not all([
            exists('{}.{}'.format(fasta_file_path, suffix))
            for suffix in ['bwt', 'pac', 'ann', 'amb', 'sa']
    ]):
        print('Could not find BWA index; building it ...')
        run_command('bwa index {}'.format(fasta_file_path))

    output_bam_file_path = fastq_gz_file_paths[0] + '.align_fastq_gzs_using_bwa.bam'

    run_command('bwa mem -t {} {} {} | samtools view -Sb -@ {} > {}'.format(
        n_jobs, fasta_file_path, ' '.join(fastq_gz_file_paths), n_jobs,
        output_bam_file_path))

    return index_bam_using_samtools(
        sort_bam_using_samtools(output_bam_file_path, n_jobs=n_jobs),
        n_jobs=n_jobs)
