from os.path import exists, join, split

from . import RESOURCE_DIRECTORY_PATH
from .print_and_run_command import print_and_run_command
from .support.support.network import download

# TODO: consider ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz


def make_reference_genome(directory_path, overwrite=False):

    final_fa_file_path = join(
        directory_path,
        'GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa')

    final_fa_gz_file_path = final_fa_file_path + '.gz'

    if not overwrite and exists(final_fa_gz_file_path):

        raise FileExistsError(final_fa_gz_file_path)

    fa_gz_file_path = join(directory_path,
                           'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz')

    download(
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/{}'.
        format(split(fa_gz_file_path)[-1]), directory_path)

    print_and_run_command('gzip --decompress --to-stdout {} > {}'.format(
        fa_gz_file_path, final_fa_file_path))

    print_and_run_command('cat {} >> {}'.format(
        join(RESOURCE_DIRECTORY_PATH, 'bwa.kit', 'resource-GRCh38',
             'hs38DH-extra.fa'), final_fa_file_path))

    print_and_run_command('gzip {} --to-stdout > {}'.format(
        final_fa_file_path, final_fa_gz_file_path))

    print_and_run_command('cp -f {} {}'.format(
        join(RESOURCE_DIRECTORY_PATH, 'bwa.kit', 'resource-GRCh38',
             'hs38DH.fa.alt'), final_fa_gz_file_path + '.alt'))

    return final_fa_gz_file_path
