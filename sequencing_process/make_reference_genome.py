from os.path import isfile, split

from . import RESOURCE_DIRECTORY_PATH
from ._print_and_run_command import _print_and_run_command
from .support.support.network import download

# TODO: consider ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz


def make_reference_genome(
        directory_path,
        overwrite=False,
):

    resource_grch_directory_path = '{}/bwa.kit/resource-GRCh38'.format(
        RESOURCE_DIRECTORY_PATH)

    final_fa_file_path = '{}/GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa'.format(
        directory_path)

    final_fa_gz_file_path = final_fa_file_path + '.gz'

    if not overwrite and isfile(final_fa_gz_file_path):

        raise FileExistsError(final_fa_gz_file_path)

    fa_gz_file_path = '{}/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz'.format(
        directory_path)

    download(
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/{}'.
        format(split(fa_gz_file_path)[-1]),
        directory_path,
    )

    _print_and_run_command('gzip --decompress --to-stdout {} > {}'.format(
        fa_gz_file_path,
        final_fa_file_path,
    ))

    _print_and_run_command('cat {}/hs38DH-extra.fa >> {}'.format(
        resource_grch_directory_path,
        final_fa_file_path,
    ))

    _print_and_run_command('gzip {} --to-stdout > {}'.format(
        final_fa_file_path,
        final_fa_gz_file_path,
    ))

    _print_and_run_command('cp -f {}/hs38DH.fa.alt {}.alt'.format(
        resource_grch_directory_path,
        final_fa_file_path,
    ))

    return final_fa_gz_file_path
