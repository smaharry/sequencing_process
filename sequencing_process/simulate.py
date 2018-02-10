from . import print_and_run_command
from .bgzip_and_tabix import bgzip_and_tabix


def simulate_sequences_using_dwgsim(fasta_file_path,
                                    output_file_path_prefix,
                                    n_sequence=1000,
                                    fraction_variant=0.001,
                                    fraction_indel_variant=0.1):
    """
    Simulate sequences using dwgsim.
    Arguments:
        fasta_file_path (str):
        output_file_path_prefix (str):
        n_sequence (int):
        fraction_variant (float):
        fraction_indel_variant (float):
    Returns:
    """

    print_and_run_command(
        'dwgsim -N {} -1 150 -2 150 -r {} -R {} {} {}'.format(
            n_sequence, fraction_variant, fraction_indel_variant,
            fasta_file_path, output_file_path_prefix))

    print_and_run_command(
        'gzip {}'.format(output_file_path_prefix + '.bwa.read1.fastq'))

    print_and_run_command(
        'gzip {}'.format(output_file_path_prefix + '.bwa.read2.fastq'))

    print_and_run_command(
        'rm -rf {}'.format(output_file_path_prefix + '.bfast.fastq'))

    bgzip_and_tabix(output_file_path_prefix + '.mutations.vcf')
