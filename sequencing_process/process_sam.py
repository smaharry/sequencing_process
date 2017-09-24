from .gz import bgzip, tabix
from .support.support.subprocess_ import run_command


def samptools_sam_to_bam(sam_file_path, n_jobs=1):
    """
    Compress .sam file to .bam file.
    Arguments:
        sam_file_path (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_bam_file_path = sam_file_path + '.samptools_sam_to_bam.bam'

    command = 'samtools view -Sb -@ {} {} > {}'.format(n_jobs, sam_file_path,
                                                       output_bam_file_path)

    run_command(command)

    return output_bam_file_path


def freebayes(bam_file_path, fasta_file_path):
    """
    Call variants on .bam file.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
    Returns:
        str:
    """

    output_vcf_file_path = bam_file_path + '.freebayes.vcf'

    command = 'freebayes --fasta-reference {} {} > {}'.format(
        fasta_file_path, bam_file_path, output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))
