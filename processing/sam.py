from .gz import bgzip, tabix


def compress_sam_to_bam(sam_file_path):
    """
    Arguments:
        sam_file_path (str):
    Returns:
        str:
    """

    output_bam_file_path = sam_file_path + '.bam'

    command =

    run_command(command)

    return output_bam_file_path


def call_variants(bam_file_path):
    """
    Arguments:
        bam_file_path (str):
    Returns:
        str:
    """

    output_vcf_file_path = bam_file_path + '.vcf'

    command =

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))
