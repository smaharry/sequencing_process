from . import CHROMOSOMES
from .bgzip_and_tabix import bgzip_and_tabix
from .process_vcf_gz import concatenate_vcf_gzs_using_bcftools
from .support.support.multiprocess import multiprocess
from .support.support.subprocess_ import run_command


def sort_bam_using_samtools(bam_file_path, n_jobs=1):
    """
    Sort .bam file using samtools.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_bam_file_path = bam_file_path.replace(
        '.bam', '.sort_bam_using_samtools.bam')

    run_command('samtools sort -@ {} {} > {}'.format(n_jobs, bam_file_path,
                                                     output_bam_file_path))

    return output_bam_file_path


def index_bam_using_samtools(bam_file_path, n_jobs=1):
    """
    Index .bam file using samtools.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
    Returns:
        str:
    """

    run_command('samtools index -@ {} {}'.format(n_jobs, bam_file_path))

    return bam_file_path


def remove_duplicates_in_bam_using_picard(bam_file_path,
                                          maximum_memory='8G',
                                          n_jobs=1):
    """
    Remove duplicates in .bam file using picard.
    Arguments:
        bam_file_path (str):
        maximum_memory (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_bam_file_path = bam_file_path.replace(
        '.bam', '.remove_duplicates_in_bam_using_picard.bam')

    run_command('picard -Xmx{} MarkDuplicates I={} O={} M={}'.format(
        maximum_memory, bam_file_path, output_bam_file_path,
        output_bam_file_path[:-4]))

    return index_bam_using_samtools(output_bam_file_path, n_jobs=n_jobs)


def call_variants_on_bam_using_freebayes_and_multiprocess(
        bam_file_path, fasta_file_path, chromosomes=CHROMOSOMES, n_jobs=2):
    """
    Call variants on .bam file using freebayes and multiprocess.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        chromosomes (iterable):
        n_jobs (int):
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
                                         n_jobs=1):
    """
    Call variants on .bam file using freebayes.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        regions (str):
        n_jobs (int):
    Returns:
        str:
    """

    additional_arguments = ''
    if regions:
        additional_arguments += '-r {}'.format(regions)

    output_vcf_gz_file_path = bam_file_path.replace(
        '.bam', '.call_variants_on_bam_using_freebayes.vcf')

    run_command('freebayes -f {} {} {} > {}'.format(
        fasta_file_path, additional_arguments, bam_file_path,
        output_vcf_gz_file_path))

    return bgzip_and_tabix(output_vcf_gz_file_path)
