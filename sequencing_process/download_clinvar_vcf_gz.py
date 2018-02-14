from os.path import exists, join

from .support.support.network import download


def download_clinvar_vcf_gz(directory_path, version=None, overwrite=False):
    """
    Download the latest clinvar.vcf.gz.
    Arguments:
        directory_path (str):
        version (str):
        overwrite (bool):
    Returns:
    """

    clinvar_vcf_gz_file_path = join(directory_path, 'clinvar.vcf.gz')
    if not overwrite and exists(clinvar_vcf_gz_file_path):
        raise FileExistsError(clinvar_vcf_gz_file_path)

    if version is None:
        version_suffix = ''
    else:
        version_suffix = '_' + version

    for url in (
            'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar{}.vcf.gz'.
            format(version_suffix),
            'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar{}.vcf.gz.tbi'.
            format(version_suffix), ):
        download(url, directory_path)
