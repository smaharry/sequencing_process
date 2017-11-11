from inspect import stack
from os import remove
from os.path import dirname, isfile, join

from . import CHROMOSOMES
from .bgzip_and_tabix import bgzip_and_tabix
from .process_vcf_gz import concatenate_vcf_gzs_using_bcftools
from .support.support.multiprocess import multiprocess
from .support.support.subprocess_ import run_command_and_monitor


def sort_bam_using_samtools(bam_file_path, n_jobs=1, overwrite=False):
    """
    Sort .bam file using samtools.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_bam_file_path = join(dirname(bam_file_path), stack()[0][3] + '.bam')

    if isfile(output_bam_file_path):
        raise FileExistsError('{} exists.'.format(output_bam_file_path))

    run_command_and_monitor(
        'samtools sort --threads {} {} > {}'.format(n_jobs, bam_file_path,
                                                    output_bam_file_path),
        print_command=True)

    remove(bam_file_path)

    return output_bam_file_path


def remove_duplicates_in_bam_using_picard(bam_file_path,
                                          maximum_memory='8G',
                                          n_jobs=1,
                                          overwrite=False):
    """
    Remove duplicates in .bam file using picard.
    Arguments:
        bam_file_path (str):
        maximum_memory (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_bam_file_path = join(dirname(bam_file_path), stack()[0][3] + '.bam')

    if isfile(output_bam_file_path):
        raise FileExistsError('{} exists.'.format(output_bam_file_path))

    metrics_file_path = output_bam_file_path + '.metrics'

    run_command_and_monitor(
        'picard -Xmx{} MarkDuplicates REMOVE_DUPLICATES=true INPUT={} OUTPUT={} METRICS_FILE={}'.
        format(maximum_memory, bam_file_path, output_bam_file_path,
               metrics_file_path),
        print_command=True)

    return index_bam_using_samtools(output_bam_file_path, n_jobs=n_jobs)


def index_bam_using_samtools(bam_file_path, n_jobs=1, overwrite=False):
    """
    Index .bam file using samtools.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    run_command_and_monitor(
        'samtools index -@ {} {}'.format(n_jobs, bam_file_path),
        print_command=True)

    return bam_file_path


def call_variants_on_bam_using_freebayes_and_multiprocess(
        bam_file_path,
        fasta_file_path,
        chromosomes=CHROMOSOMES,
        n_jobs=2,
        overwrite=False):
    """
    Call variants on .bam file using freebayes and multiprocess.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        chromosomes (iterable):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    ps = multiprocess(
        call_variants_on_bam_using_freebayes,
        [[bam_file_path, fasta_file_path, c] for c in chromosomes],
        n_jobs=n_jobs)

    return concatenate_vcf_gzs_using_bcftools(ps, n_jobs=n_jobs)


def call_variants_on_bam_using_freebayes(bam_file_path,
                                         fasta_file_path,
                                         regions=None,
                                         n_jobs=1,
                                         overwrite=False):
    """
    Call variants on .bam file using freebayes.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        regions (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_vcf_file_path = join(dirname(bam_file_path), stack()[0][3] + '.vcf')

    # if isfile(output_vcf_file_path):
    #     raise FileExistsError('{} exists.'.format(output_vcf_file_path))

    additional_argument = ''
    if regions:
        additional_argument += '--region {}'.format(regions)

    if additional_argument:
        output_vcf_file_path = output_vcf_file_path.replace(
            '.vcf', '.{}.vcf'.format(additional_argument.replace(' ', '_')))

    # run_command_and_monitor(
    #     'freebayes --fasta-reference {} {} {} > {}'.format(
    #         fasta_file_path, additional_argument, bam_file_path,
    #         output_vcf_file_path),
    #     print_command=True)

    print(bgzip_and_tabix(output_vcf_file_path, n_jobs=n_jobs))
    return bgzip_and_tabix(output_vcf_file_path, n_jobs=n_jobs)
