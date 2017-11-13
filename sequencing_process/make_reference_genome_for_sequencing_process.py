from os import chdir

from .support.support.network import download
from .support.support.subprocess_ import run_command_and_monitor


def make_reference_genome_for_sequencing_process(directory_path):
    """
    """

    full_name = 'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz'
    extra_name = 'hs38DH-extra.fa'
    full_plus_extra_name = 'GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa'

    chdir(directory_path)

    download(
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/{}'.
        format(full_name))
    download(''.format(extra_name))

    run_command_and_monitor(
        'gzip -dc {0} > {2}; cat {1} >> {2}; gzip {2}'.format(
            full_name, extra_name, full_plus_extra_name))
