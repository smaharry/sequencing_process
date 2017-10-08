from . import CHROMOSOMES
from .support.support.multiprocess import multiprocess
from .support.support.subprocess_ import run_command


def remove_duplicates_in_bam_using_picard(bam_file_path, maximum_memory='8G'):
    """
    Mark duplicates in .bam file using picard.
    Arguments:
        bam_file_path (str):
        maximum_memory (str):
    Returns:
        str:
    """

    output_bam_file_path = bam_file_path + '.remove_duplicates_in_bam_using_picard.bam'
    report_file_path = bam_file_path + '.remove_duplicates_in_bam_using_picard'

    command = 'picard -Xmx{} MarkDuplicates I={} O={} M={}'.format(
        maximum_memory, bam_file_path, output_bam_file_path, report_file_path)

    run_command(command)

    return output_bam_file_path


def call_variants_on_bam_using_multiprocess_using_freebayes(
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

    args = [[bam_file_path, fasta_file_path, '{}'.format(c)]
            for c in chromosomes]

    output_vcf_gz_file_paths = multiprocess(
        call_variants_on_bam_using_freebayes, args, n_jobs=n_jobs)

    output_vcf_gz_file_path = bam_file_path + '.call_variants_on_bam_using_multiprocess_using_freebayes.vcf.gz'

    command = 'bcftools concat -a {} --threads {} | bgzip -cf -@ {} > {}; tabix -f {}'.format(
        ' '.join(output_vcf_gz_file_paths), n_jobs, n_jobs,
        output_vcf_gz_file_path, output_vcf_gz_file_path)

    run_command(command)

    return output_vcf_gz_file_path


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

    output_vcf_gz_file_path = bam_file_path + '.call_variants_on_bam_using_freebayes.vcf.gz'

    command = 'freebayes -f {} {} {} | bgzip -cf -@ {} > {}; tabix -f {}'.format(
        fasta_file_path, additional_arguments, bam_file_path, n_jobs,
        output_vcf_gz_file_path, output_vcf_gz_file_path)

    run_command(command)

    return output_vcf_gz_file_path
