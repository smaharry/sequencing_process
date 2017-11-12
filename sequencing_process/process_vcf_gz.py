from inspect import stack
from os import remove
from os.path import dirname, exists, join

from .bgzip_and_tabix import bgzip_and_tabix
from .support.support.subprocess_ import run_command_and_monitor


def concatenate_vcf_gzs_using_bcftools(vcf_gz_file_paths,
                                       n_jobs=1,
                                       overwrite=False):
    """
    Concatenate .vcf.gz files using bcftools.
    Arguments
        vcf_gz_file_paths (iterable):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_vcf_file_path = join(
        dirname(vcf_gz_file_paths[0]), stack()[0][3] + '.vcf')

    if not overwrite and exists(output_vcf_file_path):
        raise FileExistsError(output_vcf_file_path)

    run_command_and_monitor(
        'bcftools concat --allow-overlaps --threads {} {} > {}'.format(
            n_jobs, ' '.join(vcf_gz_file_paths), output_vcf_file_path),
        print_command=True)

    for vcf_gz_file_path in vcf_gz_file_paths:
        remove(vcf_gz_file_path)
        remove(vcf_gz_file_path + '.tbi')

    return bgzip_and_tabix(
        output_vcf_file_path, n_jobs=n_jobs, overwrite=overwrite)


def extract_regions_from_vcf_gz_using_bcftools(vcf_gz_file_path,
                                               regions,
                                               n_jobs=1,
                                               overwrite=False):
    """
    Extract regions from .vcf.gz file using bcftools.
    Arguments:
        vcf_gz_file_path (str):
        regions (iterable):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_vcf_file_path = join(
        dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')

    if not overwrite and exists(output_vcf_file_path):
        raise FileExistsError(output_vcf_file_path)

    run_command_and_monitor(
        'bcftools view --regions {} --threads {} {} > {}'.format(
            ','.join(regions), n_jobs, vcf_gz_file_path, output_vcf_file_path),
        print_command=True)

    return bgzip_and_tabix(
        output_vcf_file_path, n_jobs=n_jobs, overwrite=overwrite)


def annotate_vcf_gz_using_snpeff(vcf_gz_file_path,
                                 genomic_assembly='GRCh38.86',
                                 maximum_memory='8G',
                                 n_jobs=1,
                                 overwrite=False):
    """
    Annotate .vcf.gz file using snpeff.
    Arguments:
        vcf_gz_file_path (str):
        genomic_assembly (str):
        maximum_memory (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_vcf_file_path = join(
        dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')

    if not overwrite and exists(output_vcf_file_path):
        raise FileExistsError(output_vcf_file_path)

    run_command_and_monitor(
        'snpEff -Xmx{} -htmlStats {}.stats.html -csvStats {}.stats.csv -verbose -noLog {} {} > {}'.
        format(maximum_memory, output_vcf_file_path, output_vcf_file_path,
               genomic_assembly, vcf_gz_file_path, output_vcf_file_path),
        print_command=True)

    return bgzip_and_tabix(
        output_vcf_file_path, n_jobs=n_jobs, overwrite=overwrite)


def annotate_vcf_gz_using_bcftools(vcf_gz_file_path,
                                   annotation_file_path,
                                   additional_arguments=['--columns =ID,INFO'],
                                   n_jobs=1,
                                   overwrite=False):
    """
    Annotate .vcf.gz file using bcftools.
    Arguments:
        vcf_gz_file_path (str):
        annotation_file_path (str):
        additional_arguments (list):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_vcf_file_path = join(
        dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')

    if not overwrite and exists(output_vcf_file_path):
        raise FileExistsError(output_vcf_file_path)

    run_command_and_monitor(
        'bcftools annotate --annotations {} --threads {} {} {} > {}'.format(
            annotation_file_path, n_jobs, ' '.join(additional_arguments),
            vcf_gz_file_path, output_vcf_file_path),
        print_command=True)

    return bgzip_and_tabix(
        output_vcf_file_path, n_jobs=n_jobs, overwrite=overwrite)


def filter_vcf_gz_using_bcftools(
        vcf_gz_file_path,
        include_expression='10 <= DP & 10 <= QUAL & 10 <= QUAL / AO & 1 <= SRF & 1 <= SRR & 1 <= SAF & 1 <= SAR & 1 <= RPR & 1 <= RPL',
        n_jobs=1,
        overwrite=False):
    """
    Filter .vcf.gz file using bcftools.
    Arguments:
        vcf_gz_file_path (str):
        include_expression (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_vcf_file_path = join(
        dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')

    if not overwrite and exists(output_vcf_file_path):
        raise FileExistsError(output_vcf_file_path)

    run_command_and_monitor(
        'bcftools view --include \'{}\' --threads {} {} > {}'.format(
            include_expression, n_jobs, vcf_gz_file_path,
            output_vcf_file_path),
        print_command=True)

    return bgzip_and_tabix(
        output_vcf_file_path, n_jobs=n_jobs, overwrite=overwrite)
