from os.path import isfile

from .support.support.network import download


def download_clinvar_vcf_gz(
        directory_path,
        version=None,
        overwrite=False,
):

    if version is None:

        clinvar_vcf_gz_file_name = 'clinvar.vcf.gz'

    else:

        clinvar_vcf_gz_file_name = 'clinvar_{}.vcf.gz'.format(version)

    clinvar_vcf_gz_file_path = '{}/{}'.format(
        directory_path,
        clinvar_vcf_gz_file_name,
    )

    if not overwrite and isfile(clinvar_vcf_gz_file_path):

        raise FileExistsError(clinvar_vcf_gz_file_path)

    for url in (
            'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/{}'.format(
                clinvar_vcf_gz_file_name),
            'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/{}.tbi'.format(
                clinvar_vcf_gz_file_name),
    ):

        download(
            url,
            directory_path,
        )

    return clinvar_vcf_gz_file_path
