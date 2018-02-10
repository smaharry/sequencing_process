from inspect import stack
from os.path import dirname, exists, join

from . import print_and_run_command
from .bgzip_and_tabix import bgzip_and_tabix
from .process_vcf_gz import concatenate_vcf_gzs_using_bcftools_concat
from .support.support.multiprocess import multiprocess


def sort_and_index_bam_using_samtools_sort_and_index(
        bam_file_path,
        remove_input_bam_file_path=False,
        n_job=1,
        output_bam_file_path=None,
        overwrite=False):
    """
    Sort and index .bam file using samtools sort and index.
    Arguments:
        bam_file_path (str):
        remove_input_bam_file_path (bool):
        n_job (int):
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

    print_and_run_command('samtools sort --threads {} {} > {}'.format(
        n_job, bam_file_path, output_bam_file_path))

    if remove_input_bam_file_path:
        print_and_run_command('rm -rf {}'.format(bam_file_path))

    return index_bam_using_samtools_index(
        output_bam_file_path, n_job=n_job, overwrite=overwrite)


def index_bam_using_samtools_index(bam_file_path, n_job=1, overwrite=False):
    """
    Index .bam file using samtools.
    Arguments:
        bam_file_path (str):
        n_job (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_bai_file_path = bam_file_path + '.bai'
    if not overwrite and exists(output_bai_file_path):
        raise FileExistsError(output_bai_file_path)

    print_and_run_command(
        'samtools index -@ {} {}'.format(n_job, bam_file_path))

    return bam_file_path


def mark_duplicates_in_bam_using_picard_markduplicates(
        bam_file_path,
        maximum_memory='8G',
        remove_duplicates=False,
        remove_input_bam_file_path_and_its_index=False,
        n_job=1,
        output_bam_file_path=None,
        overwrite=False):
    """
    Remove duplicates in .bam file using picard markduplicates.
    Arguments:
        bam_file_path (str):
        maximum_memory (str):
        remove_duplicates (bool):
        remove_input_bam_file_path_and_its_index (bool):
        n_job (int):
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

    metrics_file_path = output_bam_file_path + '.metrics'
    print_and_run_command(
        'picard -Xmx{} MarkDuplicates REMOVE_DUPLICATES={} INPUT={} OUTPUT={} METRICS_FILE={}'.
        format(maximum_memory,
               str(remove_duplicates).lower(), bam_file_path,
               output_bam_file_path, metrics_file_path))

    if remove_input_bam_file_path_and_its_index:
        print_and_run_command('rm -rf {}'.format(bam_file_path))
        print_and_run_command('rm -rf {}'.format(bam_file_path + '.bai'))

    print('{}:'.format(metrics_file_path))
    with open(metrics_file_path) as file_:
        print(file_.read())

    return index_bam_using_samtools_index(
        output_bam_file_path, n_job=n_job, overwrite=overwrite)


def check_bam_using_samtools_flagstat(bam_file_path,
                                      n_job=1,
                                      output_file_path=None,
                                      overwrite=False):
    """
    Check .bam file using samtools flagstat.
    Arguments:
        bam_file_path (str):
        n_job (int):
        output_file_path (str):
        overwrite (bool):
    Returns:
    """

    if not output_file_path:
        output_file_path = bam_file_path + '.flagstat'
    if not overwrite and exists(output_file_path):
        raise FileExistsError(output_file_path)

    print_and_run_command('samtools flagstat --threads {} {} > {}'.format(
        n_job, bam_file_path, output_file_path))

    print('{}:'.format(output_file_path))
    with open(output_file_path) as file_:
        print(file_.read())


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
    """

    plot_zip_prefix_path = fastq_gz_or_bam_file_path + '.plot'
    plot_tsv_file_path = plot_zip_prefix_path + '.tsv'
    if not overwrite:
        if exists(plot_zip_prefix_path + '.zip'):
            raise FileExistsError(plot_zip_prefix_path + '.zip')
        if exists(plot_tsv_file_path):
            raise FileExistsError(plot_tsv_file_path)

    print_and_run_command('fastqp --kmer {} --output {} --text {} {}'.format(
        kmer_length, plot_zip_prefix_path, plot_tsv_file_path,
        fastq_gz_or_bam_file_path))


def call_variants_on_bam_using_freebayes_and_multiprocess(
        bam_file_path,
        fasta_file_path,
        regions,
        n_job=2,
        output_vcf_file_path=None,
        overwrite=False):
    """
    Call variants on .bam file using freebayes via multiprocess.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        regions (iterable):
        n_job (int):
        output_vcf_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    output_vcf_gz_file_path = concatenate_vcf_gzs_using_bcftools_concat(
        multiprocess(
            call_variants_on_bam_using_freebayes, ((
                bam_file_path,
                fasta_file_path,
                region,
                1,
                None,
                overwrite, ) for region in regions),
            n_job=n_job),
        remove_input_vcf_gz_file_paths_and_their_indices=True,
        n_job=n_job,
        output_vcf_file_path=output_vcf_file_path,
        overwrite=overwrite)

    return output_vcf_gz_file_path


def call_variants_on_bam_using_freebayes(bam_file_path,
                                         fasta_file_path,
                                         regions=None,
                                         n_job=1,
                                         output_vcf_file_path=None,
                                         overwrite=False):
    """
    Call variants on .bam file using freebayes.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        regions (str):
        n_job (int):
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

    print_and_run_command('freebayes --fasta-reference {} {} {} > {}'.format(
        fasta_file_path, ' '.join(additional_arguments), bam_file_path,
        output_vcf_file_path))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)
