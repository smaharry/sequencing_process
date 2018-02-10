from os.path import exists, join

from .support.support.network import download


def download_clinvar_vcf_gz(directory_path, overwrite=False):
    """
    Download the latest clinvar.vcf.gz.
    Arguments:
        directory_path (str):
        overwrite (bool):
    Returns:
    """

    clinvar_vcf_gz_file_path = join(directory_path, 'clinvar.vcf.gz')
    if not overwrite and exists(clinvar_vcf_gz_file_path):
        raise FileExistsError(clinvar_vcf_gz_file_path)

    for url in (
            'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz',
            'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi',
    ):
        download(url, directory_path)
