from inspect import stack
from os import remove
from os.path import dirname, exists, join

from .bgzip_and_tabix import bgzip_and_tabix
from .process_vcf_gz import concatenate_vcf_gzs_using_bcftools
from .support.support.multiprocess import multiprocess
from .support.support.subprocess_ import run_command


def sort_and_index_bam_using_samtools(bam_file_path,
                                      n_jobs=1,
                                      output_bam_file_path=None,
                                      overwrite=False):
    """
    Sort and index .bam file using samtools.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
        output_bam_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_bam_file_path:
        output_bam_file_path = join(
            dirname(bam_file_path), stack()[0][3] + '.bam')

    if not overwrite and exists(output_bam_file_path):
        raise FileExistsError(output_bam_file_path)

    run_command(
        'samtools sort --threads {} {} > {}'.format(n_jobs, bam_file_path,
                                                    output_bam_file_path),
        print_command=True)

    print('Consider removing unsorted .bam file {}.'.format(bam_file_path))

    return index_bam_using_samtools_index(
        output_bam_file_path, n_jobs=n_jobs, overwrite=overwrite)


def index_bam_using_samtools_index(bam_file_path, n_jobs=1, overwrite=False):
    """
    Index .bam file using samtools.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_bai_file_path = bam_file_path + '.bai'

    if not overwrite and exists(output_bai_file_path):
        raise FileExistsError(output_bai_file_path)

    run_command(
        'samtools index -@ {} {}'.format(n_jobs, bam_file_path),
        print_command=True)

    return bam_file_path


def mark_duplicates_in_bam_using_picard(bam_file_path,
                                        maximum_memory='8G',
                                        remove=True,
                                        n_jobs=1,
                                        output_bam_file_path=None,
                                        overwrite=False):
    """
    Remove duplicates in .bam file using picard.
    Arguments:
        bam_file_path (str):
        maximum_memory (str):
        remove (bool):
        n_jobs (int):
        output_bam_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_bam_file_path:
        output_bam_file_path = join(
            dirname(bam_file_path), stack()[0][3] + '.bam')

    if not overwrite and exists(output_bam_file_path):
        raise FileExistsError(output_bam_file_path)

    run_command(
        'picard -Xmx{} MarkDuplicates REMOVE_DUPLICATES={} INPUT={} OUTPUT={} METRICS_FILE={}.metrics'.
        format(maximum_memory,
               str(remove).lower(), bam_file_path, output_bam_file_path,
               output_bam_file_path),
        print_command=True)

    return index_bam_using_samtools_index(
        output_bam_file_path, n_jobs=n_jobs, overwrite=overwrite)


def check_fastq_gz_or_bam_using_fastqp(fastq_gz_or_bam_file_path,
                                       kmer_length=7,
                                       overwrite=False):
    """
    Check .fastq.gz or .bam file using fastqp.
    Arguments:
        fastq_gz_or_bam_file_path (str):
        kmer_length (int):
        overwrite (bool):
    Returns:
        None
    """

    plot_zip_prefix_path = fastq_gz_or_bam_file_path + '.plot'
    plot_tsv_file_path = plot_zip_prefix_path + '.tsv'

    if not overwrite:
        if exists(plot_zip_prefix_path + '.zip'):
            raise FileExistsError(plot_zip_prefix_path + '.zip')
        if exists(plot_tsv_file_path):
            raise FileExistsError(plot_tsv_file_path)

    run_command(
        'fastqp --kmer {} --output {} --text {} --count-duplicates True {}'.
        format(kmer_length, plot_zip_prefix_path, plot_tsv_file_path,
               fastq_gz_or_bam_file_path),
        print_command=True)


def check_bam_using_samtools_flagstat(bam_file_path,
                                      n_jobs=1,
                                      output_file_path=None,
                                      overwrite=False):
    """
    Check .bam file using samtools flagstat.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
        output_file_path (str):
        overwrite (bool):
    Returns:
        None
    """

    if not output_file_path:
        output_file_path = join(dirname(bam_file_path), stack()[0][3])

    if not overwrite and exists(output_file_path):
        raise FileExistsError(output_file_path)

    run_command(
        'samtools flagstat --threads {} {} > {}'.format(
            n_jobs, bam_file_path, output_file_path),
        print_command=True)


def call_variants_on_bam_using_freebayes_and_multiprocess(
        bam_file_path,
        fasta_file_path,
        regions,
        n_jobs=2,
        output_vcf_file_path=None,
        overwrite=False):
    """
    Call variants on .bam file using freebayes and multiprocess.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        regions (iterable):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    ps = multiprocess(
        call_variants_on_bam_using_freebayes,
        [[bam_file_path, fasta_file_path, r, 1, None, overwrite]
         for r in regions],
        n_jobs=n_jobs)

    output_vcf_gz_file_path = concatenate_vcf_gzs_using_bcftools(
        ps,
        n_jobs=n_jobs,
        output_vcf_file_path=output_vcf_file_path,
        overwrite=overwrite)

    for p in ps:
        remove(p)
        remove(p + '.tbi')

    return output_vcf_gz_file_path


def call_variants_on_bam_using_freebayes(bam_file_path,
                                         fasta_file_path,
                                         regions=None,
                                         n_jobs=1,
                                         output_vcf_file_path=None,
                                         overwrite=False):
    """
    Call variants on .bam file using freebayes.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        regions (str):
        n_jobs (int):
        output_vcf_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_vcf_file_path:
        output_vcf_file_path = join(
            dirname(bam_file_path), stack()[0][3] + '.vcf')

    additional_arguments = []

    if regions:
        additional_arguments.append('--region {}'.format(regions))

    if any(additional_arguments):
        output_vcf_file_path = output_vcf_file_path.replace(
            '.vcf',
            '.{}.vcf'.format(' '.join(additional_arguments).replace(' ', '_')))

    if not overwrite and exists(output_vcf_file_path):
        raise FileExistsError(output_vcf_file_path)

    run_command(
        'freebayes --fasta-reference {} {} {} > {}'.format(
            fasta_file_path, ' '.join(additional_arguments), bam_file_path,
            output_vcf_file_path),
        print_command=True)

    return bgzip_and_tabix(
        output_vcf_file_path, n_jobs=n_jobs, overwrite=overwrite)
