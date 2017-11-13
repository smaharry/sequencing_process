from os import chdir

from .support.support.network import download
from .support.support.subprocess_ import run_command_and_monitor


def make_reference_genome_for_sequencing_process(directory_path):
    """
    Make make reference genome for sequencing process.
    Get NCBI's GRCH38 sequences designed for sequencing process (0).
    Get Heng Li's decoy and HLA sequences (1).
    Concatenate 0 and 1 into 2 and gzip 2 together.
    Copy 1.alt to 2.alt because bwa needs 2.alt to mark a aliignment as alt and
    align alt sequences to non-alt sequences.
    Arguments:
        directory_path (str):
    Returns:
        None
    """

    chdir(directory_path)

    full_name = 'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz'
    extra_name = 'hs38DH-extra.fa'
    full_plus_extra_name = 'GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa'

    download(
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/{}'.
        format(full_name))

    download(
        'https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download'
    )
    run_command_and_monitor('gzip -dc download | tar -xf -')

    run_command_and_monitor(
        'gzip -dc {0} > {2}; cat {1} >> {2}; gzip {2}'.format(
            full_name, extra_name, full_plus_extra_name))


if __name__ == '__main__':
    make_reference_genome_for_sequencing_process()
