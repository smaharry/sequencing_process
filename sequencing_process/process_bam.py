from . import CHROMOSOMES
from .support.support.multiprocess import multiprocess
from .support.support.subprocess_ import run_command


def mark_duplicates_in_bam_with_samtools(bam_file_path, n_jobs):
    """
    Mark duplicates in .bam file with samtools.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_bam_file_path = bam_file_path + '.mark_duplicates_in_bam_with_samtools.bam'

    # command = 'samtools sort -n -@ {} {} | samtools fixmate -m -@ {} | samtools sort -@ {} | samtools markdup -@ {} > {}; samtools index -@ {} {}'.format(
    #     n_jobs, bam_file_path, n_jobs, n_jobs, n_jobs, output_bam_file_path,
    #     n_jobs, output_bam_file_path)
    command = 'samtools rmdup -@ {} {} > {}'.format(n_jobs, bam_file_path,
                                                    output_bam_file_path)

    run_command(command)

    return output_bam_file_path


def call_variants_on_bam_using_multiprocess_with_freebayes(
        bam_file_path, fasta_file_path, n_jobs, chromosomes=CHROMOSOMES):
    """
    Call variants on .bam file using multiprocess with freebayes.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        n_jobs (int):
        chromosomes (iterable):
    Returns:
        str:
    """

    args = [[bam_file_path, fasta_file_path, '{}'.format(c)]
            for c in chromosomes]

    returns = multiprocess(
        call_variants_on_bam_with_freebayes, args, n_jobs=n_jobs)

    output_vcf_gz_file_path = bam_file_path + '.call_variants_on_bam_using_multiprocess_with_freebayes.vcf.gz'

    command = 'bcftools concat {} | bgzip -fc -@ {} > {}; tabix -f {}'.format(
        ' '.join(returns), n_jobs, output_vcf_gz_file_path,
        output_vcf_gz_file_path)

    run_command(command)

    return output_vcf_gz_file_path


def call_variants_on_bam_with_freebayes(bam_file_path,
                                        fasta_file_path,
                                        regions=None,
                                        n_jobs=1):
    """
    Call variants on .bam file with freebayes.
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

    if additional_arguments:
        output_vcf_gz_file_path = bam_file_path + '.call_variants_on_bam_with_freebayes.{}.vcf.gz'.format(
            additional_arguments.replace(' ', '_'))

    else:
        output_vcf_gz_file_path = bam_file_path + '.call_variants_on_bam_with_freebayes.vcf.gz'

    command = 'freebayes -f {} {} {} | bgzip -fc -@ {} > {}; tabix -f {}'.format(
        fasta_file_path, additional_arguments, bam_file_path, n_jobs,
        output_vcf_gz_file_path, output_vcf_gz_file_path)

    run_command(command)

    return output_vcf_gz_file_path
