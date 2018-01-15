from os.path import exists, join, split

from . import RESOURCE_DIRECTORY_PATH
from .support.support.network import download
from .support.support.subprocess_ import run_command


def make_reference_genome(directory_path, overwrite=False):
    """
    Make reference genome.
    Get NCBI's GRCh38 .fasta.gz file designed for sequencing process (0).
    Get Heng Li's .fasta file containing decoy and HLA sequences (1).
    Concatenate 0 and 1 into 2.
    Copy 1.alt to 2.alt (because bwa needs 2.alt to mark an alignment as alt
    and align alt sequences to non-alt sequences.)
    Arguments:
        directory_path (str):
        overwrite (bool):
    Returns:
    """

    f_e_fa_gz_file_path = join(
        directory_path,
        'GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa.gz')
    if not overwrite and exists(f_e_fa_gz_file_path):
        raise FileExistsError(f_e_fa_gz_file_path)

    f_fa_gz_file_path = join(
        directory_path, 'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz')
    download(
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/{}'.
        format(split(f_fa_gz_file_path)[-1]), directory_path)

    e_fa_file_path = join(RESOURCE_DIRECTORY_PATH, 'bwa.kit',
                          'resource-GRCh38', 'hs38DH-extra.fa')
    run_command(
        'gzip -dc {0} > {2} && cat {1} >> {2} && gzip {2}'.format(
            f_fa_gz_file_path, e_fa_file_path, f_e_fa_gz_file_path[:-3]),
        print_command=True)

    fa_alt_file_path = join(RESOURCE_DIRECTORY_PATH, 'bwa.kit',
                            'resource-GRCh38', 'hs38DH.fa.alt')
    run_command(
        'cp -f {} {}'.format(fa_alt_file_path, f_e_fa_gz_file_path + '.alt'),
        print_command=True)
