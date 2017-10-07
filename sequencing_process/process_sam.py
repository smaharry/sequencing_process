from .support.support.subprocess_ import run_command


def call_variants_on_bam(bam_file_path, fasta_file_path, n_jobs=1):
    """
    Call variants on .bam file with freebayes.
    Arguments:
        bam_file_path (str):
        fasta_file_path (str): reference .fasta file
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_file_path = bam_file_path + '.freebayes.vcf'

    # command = 'freebayes --fasta-reference {} {} | bgzip -fc -@ {} > {}; tabix -f {}'.format(
    command = 'freebayes -f {} {} > {}; bgzip -f -@ {}; tabix -f {}'.format(
        fasta_file_path,
        bam_file_path,
        output_vcf_file_path,
        n_jobs,
        output_vcf_file_path + '.gz',
        output_vcf_file_path + '.gz', )

    run_command(command)

    return output_vcf_file_path + '.gz'
