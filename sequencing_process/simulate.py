from . import print_and_run_command
from .bgzip_and_tabix import bgzip_and_tabix


def simulate_sequences_using_dwgsim(fasta_file_path,
                                    file_path_prefix,
                                    n_sequences=1000,
                                    fraction_variant=0.001,
                                    fraction_indel_variant=0.1):
    """
    Simulate sequences using dwgsim.
    Arguments:
        fasta_file_path (str):
        file_path_prefix (str):
        n_sequences (int):
        n_sequences (int):
        fraction_variant (float):
        fraction_indel_variant (float):
    Returns:
    """

    print_and_run_command(
        'dwgsim -N {} -1 150 -2 150 -r {} -R {} {} {}'.format(
            n_sequences, fraction_variant, fraction_indel_variant,
            fasta_file_path, file_path_prefix))

    print_and_run_command(
        'gzip {}'.format(file_path_prefix + '.bwa.read1.fastq'))

    print_and_run_command(
        'gzip {}'.format(file_path_prefix + '.bwa.read2.fastq'))

    print_and_run_command(
        'rm -rf {}'.format(file_path_prefix + '.bfast.fastq'))

    bgzip_and_tabix(file_path_prefix + '.mutations.vcf')
