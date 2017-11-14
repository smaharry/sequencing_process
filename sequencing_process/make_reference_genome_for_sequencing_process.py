from os import chdir
from os.path import exists, expanduser
from shutil import copy

from .support.support.network import download
from .support.support.subprocess_ import run_command_and_monitor


def make_reference_genome_for_sequencing_process(
        directory_path=expanduser('~/Downloads')):
    """
    Make reference genome for sequencing process.
    Get NCBI's GRCh38 sequences designed for sequencing process (0).
    Get Heng Li's decoy and HLA sequences (1).
    Concatenate 0 and 1 into 2 and gzip it.
    Copy 1.alt to 2.alt because bwa needs 2.alt to mark a aliignment as alt and
    align alt sequences to non-alt sequences.
    Arguments:
        directory_path (str):
    Returns:
        None
    """

    chdir(directory_path)

    full_file_name = 'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz'
    if not exists(full_file_name):
        print('Downloading {} ...'.format(full_file_name))
        download(
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/{}'.
            format(full_file_name))

    print('Downloading and uncompressing bwa.kit ...')
    download(
        'https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download'
    )
    run_command_and_monitor('gzip -dc download | tar -xf -; rm -rf download')

    extra_file_path = 'bwa.kit/resource-GRCh38/hs38DH-extra.fa'
    full_plus_extra_file_name = 'GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa'
    print('Concatenating {} and {} into {} and gzipping it ...'.format(
        full_file_name, extra_file_path, full_plus_extra_file_name))
    run_command_and_monitor(
        'gzip -dc {0} > {2}; cat {1} >> {2}; gzip {2}'.format(
            full_file_name, extra_file_path, full_plus_extra_file_name))

    alt_file_path = 'bwa.kit/resource-GRCh38/hs38DH.fa.alt'
    print('Copying {} to {}.gz.alt ...'.format(alt_file_path,
                                               full_plus_extra_file_name))
    copy(alt_file_path, '{}.alt'.format(full_plus_extra_file_name))
