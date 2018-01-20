from .bgzip_and_tabix import bgzip_and_tabix
from .support.support.path import establish_path
from .support.support.subprocess_ import run_command


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

    establish_path(file_path_prefix, 'file')

    run_command(
        'dwgsim -N {} -1 150 -2 150 -r {} -R {} {} {}'.format(
            n_sequences, fraction_variant, fraction_indel_variant,
            fasta_file_path, file_path_prefix),
        print_command=True)

    run_command(
        'gzip {}'.format(file_path_prefix + '.bwa.read1.fastq'),
        print_command=True)

    run_command(
        'gzip {}'.format(file_path_prefix + '.bwa.read2.fastq'),
        print_command=True)

    run_command(
        'rm -rf {}'.format(file_path_prefix + '.bfast.fastq'),
        print_command=True)

    bgzip_and_tabix(file_path_prefix + '.mutations.vcf')


# TODO: evaluate processing of the simulated sequences
# rows = []
# with open('1k/dwgsim_eval.tsv') as file_:
#     for line in file_:
#         line = line.strip()
#
#         if line.startswith('#'):
#             print(line)
#         else:
#             rows.append(([float(x.strip()) for x in line.split() if x]))
#
# df = DataFrame(
#     rows[:-1],
#     columns=(
#         'thr',
#         'mc',
#         'mi',
#         'mu',
#         'um',
#         'uu',
#         'n_mapped',
#         'mc^',
#         'mi^',
#         'mu^',
#         'um^',
#         'uu^',
#         'n_mapped^',
#         's',
#         'ppv',
#         'fdr',
#         's^',
#         'ppv^',
#         'fdr^', )).set_index('thr')
#
# return df
