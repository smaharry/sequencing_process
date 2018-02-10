from os.path import exists, join, split

from . import RESOURCE_DIRECTORY_PATH, print_and_run_command
from .support.support.network import download


def make_reference_genome(directory_path, overwrite=False):
    """
    Make reference genome.
    Get NCBI's GRCh38 .fasta.gz file designed for sequencing process (0).
    Get Heng Li's .fasta file containing decoy and HLA sequences (1).
    Concatenate 0 and 1 into 2.
    Copy 1.alt to 2.alt (because bwa needs 2.alt to mark an alignment as alt).
    Arguments:
        directory_path (str):
        overwrite (bool):
    Returns:
    """

    final_fa_gz_file_path = join(
        directory_path,
        'GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa.gz')
    if not overwrite and exists(final_fa_gz_file_path):
        raise FileExistsError(final_fa_gz_file_path)

    fa_gz_file_path = join(directory_path,
                           'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz')
    download(
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/{}'.
        format(split(fa_gz_file_path)[-1]), directory_path)

    print_and_run_command(
        'gzip -dc {0} > {2} && cat {1} >> {2} && gzip {2}'.format(
            fa_gz_file_path,
            join(RESOURCE_DIRECTORY_PATH, 'bwa.kit', 'resource-GRCh38',
                 'hs38DH-extra.fa'), final_fa_gz_file_path[:-3]))

    print_and_run_command('cp -f {} {}'.format(
        join(RESOURCE_DIRECTORY_PATH, 'bwa.kit', 'resource-GRCh38',
             'hs38DH.fa.alt'), final_fa_gz_file_path + '.alt'))
