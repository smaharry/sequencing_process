from inspect import stack
from os.path import dirname, exists, join

from . import RESOURCE_DIRECTORY_PATH, print_and_run_command
from .bgzip_and_tabix import bgzip_and_tabix


def concatenate_vcf_gzs_using_bcftools_concat(
        vcf_gz_file_paths,
        remove_input_vcf_gz_file_paths_and_their_indices=False,
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False):
    """
    Concatenate .vcf.gz files using bcftools concat.
    Arguments
        vcf_gz_file_paths (iterable):
        remove_input_vcf_gz_file_paths_and_their_indices (bool):
        n_job (int):
        output_vcf_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_vcf_file_path:
        output_vcf_file_path = join(
            dirname(vcf_gz_file_paths[0]), stack()[0][3] + '.vcf')
    if not overwrite and exists(output_vcf_file_path + '.gz'):
        raise FileExistsError(output_vcf_file_path + '.gz')

    print_and_run_command(
        'bcftools concat --allow-overlaps --threads {} {} > {}'.format(
            n_job, ' '.join(vcf_gz_file_paths), output_vcf_file_path))

    if remove_input_vcf_gz_file_paths_and_their_indices:
        for vcf_gz_file_path in vcf_gz_file_paths:
            print_and_run_command('rm -rf {}'.format(vcf_gz_file_path))
            print_and_run_command(
                'rm -rf {}'.format(vcf_gz_file_path + '.tbi'))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)


def rename_chromosomes_of_vcf_gz_using_bcftools_annotate(
        vcf_gz_file_path,
        map_file_path=join(RESOURCE_DIRECTORY_PATH, 'chrn_n.tsv'),
        remove_input_vcf_gz_file_path_and_its_index=False,
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False):
    """
    Rename chromosomes of .vcf.gz file using bcftools annotate.
    Arguments:
        vcf_gz_file_path (str):
        map_file_path (str):
        remove_input_vcf_gz_file_path_and_its_index (bool):
        n_job (int):
        output_vcf_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_vcf_file_path:
        output_vcf_file_path = join(
            dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')
    if not overwrite and exists(output_vcf_file_path + '.gz'):
        raise FileExistsError(output_vcf_file_path + '.gz')

    print_and_run_command(
        'bcftools annotate --rename-chrs {} --threads {} {} > {}'.format(
            map_file_path, n_job, vcf_gz_file_path, output_vcf_file_path))

    if remove_input_vcf_gz_file_path_and_its_index:
        print_and_run_command('rm -rf {}'.format(vcf_gz_file_path))
        print_and_run_command('rm -rf {}'.format(vcf_gz_file_path + '.tbi'))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)


def extract_regions_from_vcf_gz_using_bcftools_view(vcf_gz_file_path,
                                                    regions,
                                                    n_job=1,
                                                    output_vcf_file_path=None,
                                                    overwrite=False):
    """
    Extract regions from .vcf.gz file using bcftools view.
    Arguments:
        vcf_gz_file_path (str):
        regions (iterable):
        n_job (int):
        output_vcf_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_vcf_file_path:
        output_vcf_file_path = join(
            dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')

    if not overwrite and exists(output_vcf_file_path + '.gz'):
        raise FileExistsError(output_vcf_file_path + '.gz')

    print_and_run_command(
        'bcftools view --regions {} --threads {} {} > {}'.format(
            ','.join(regions), n_job, vcf_gz_file_path, output_vcf_file_path))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)


def annotate_vcf_gz_using_snpeff(
        vcf_gz_file_path,
        genomic_assembly,
        maximum_memory='12G',
        remove_input_vcf_gz_file_path_and_its_index=False,
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False):
    """
    Annotate .vcf.gz file using snpeff.
    Arguments:
        vcf_gz_file_path (str):
        genomic_assembly (str): GRCh38.86 | GRCh37.75 | hg19
        maximum_memory (str):
        remove_input_vcf_gz_file_path_and_its_index (bool):
        n_job (int):
        output_vcf_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_vcf_file_path:
        output_vcf_file_path = join(
            dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')
    if not overwrite and exists(output_vcf_file_path + '.gz'):
        raise FileExistsError(output_vcf_file_path + '.gz')

    print_and_run_command(
        'snpEff -Xmx{} -htmlStats {}.stats.html -csvStats {}.stats.csv -t -verbose -noLog {} {} > {}'.
        format(maximum_memory, output_vcf_file_path, output_vcf_file_path,
               genomic_assembly, vcf_gz_file_path, output_vcf_file_path))

    if remove_input_vcf_gz_file_path_and_its_index:
        print_and_run_command('rm -rf {}'.format(vcf_gz_file_path))
        print_and_run_command('rm -rf {}'.format(vcf_gz_file_path + '.tbi'))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)


def annotate_vcf_gz_using_bcftools_annotate(
        vcf_gz_file_path,
        annotation_file_path,
        additional_arguments=('--columns =ID,INFO', ),
        remove_input_vcf_gz_file_path_and_its_index=False,
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False):
    """
    Annotate .vcf.gz file using bcftools annotate.
    Arguments:
        vcf_gz_file_path (str):
        annotation_file_path (str):
        additional_arguments (list):
        remove_input_vcf_gz_file_path_and_its_index (bool):
        n_job (int):
        output_vcf_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_vcf_file_path:
        output_vcf_file_path = join(
            dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')
    if not overwrite and exists(output_vcf_file_path + '.gz'):
        raise FileExistsError(output_vcf_file_path + '.gz')

    print_and_run_command(
        'bcftools annotate --annotations {} --threads {} {} {} > {}'.format(
            annotation_file_path, n_job, ' '.join(additional_arguments),
            vcf_gz_file_path, output_vcf_file_path))

    if remove_input_vcf_gz_file_path_and_its_index:
        print_and_run_command('rm -rf {}'.format(vcf_gz_file_path))
        print_and_run_command('rm -rf {}'.format(vcf_gz_file_path + '.tbi'))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)


def filter_vcf_gz_using_bcftools_view(
        vcf_gz_file_path,
        include_expression='10<DP & 30<QUAL & 10<(QUAL/AO) & 1<SRF & 1<SRR & 1<SAF & 1<SAR & 1<RPR & 1<RPL',
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False):
    """
    Filter .vcf.gz file using bcftools annotate.
    Arguments:
        vcf_gz_file_path (str):
        include_expression (str):
        n_job (int):
        output_vcf_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    if not output_vcf_file_path:
        output_vcf_file_path = join(
            dirname(vcf_gz_file_path), stack()[0][3] + '.vcf')
    if not overwrite and exists(output_vcf_file_path + '.gz'):
        raise FileExistsError(output_vcf_file_path + '.gz')

    print_and_run_command(
        'bcftools view --include \'{}\' --threads {} {} > {}'.format(
            include_expression, n_job, vcf_gz_file_path, output_vcf_file_path))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)
