from . import CHROMOSOMES
from .support.support.multiprocess import multiprocess
from .support.support.subprocess_ import run_command


def call_variants_on_bam_with_multiprocess(bam_file_path,
                                           fasta_file_path,
                                           n_jobs,
                                           chromosomes=CHROMOSOMES):
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

    args = [[bam_file_path, fasta_file_path, '-r {}'.format(c)]
            for c in chromosomes]

    multiprocess(call_variants_on_bam, args, n_jobs=n_jobs)

    return


def call_variants_on_bam(bam_file_path,
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

    output_vcf_gz_file_path = bam_file_path + '.freebayes.vcf.gz'

    additional_arguments = ''

    if regions:
        additional_arguments += '-r {}'.format(regions)

    command = 'freebayes --fasta-reference {} {} {} | bgzip -fc -@ {} > {}; tabix -f {}'.format(
        fasta_file_path,
        additional_arguments,
        bam_file_path,
        n_jobs,
        output_vcf_gz_file_path,
        output_vcf_gz_file_path, )

    run_command(command)

    return output_vcf_gz_file_path
