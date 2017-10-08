from os.path import join, split

from . import CHROMOSOMES
from .process_gz import bgzip, tabix
from .support.support.subprocess_ import run_command


def picard_liftovervcf(vcf_file_path, assembly_chain_file_path,
                       target_assembly_file_path):
    """
    Re-map .vcf file coordinates.
    Arguments:
        vcf_file_path (str):
        assembly_chain_file_path (str):
        target_assembly_file_path (str):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.picard_liftovervcf.vcf'

    rejected_vcf_file_path = vcf_file_path + '.picard_liftovervcf.rejected.vcf'

    command = 'picard LiftoverVcf INPUT={} CHAIN={} REFERENCE_SEQUENCE={} OUTPUT={} REJECT={}'.format(
        vcf_file_path, assembly_chain_file_path, target_assembly_file_path,
        output_vcf_file_path, rejected_vcf_file_path)

    run_command(command)

    tabix(bgzip(rejected_vcf_file_path))

    return tabix(bgzip(output_vcf_file_path))


def bcftools_isec(vcf_file_path_1, vcf_file_path_2, n_jobs=1):
    """
    Get the intersection between 2 .vcf files.
    Arguments:
        vcf_file_path_1 (str):
        vcf_file_path_2 (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_directory_path = join(split(vcf_file_path_1)[0], 'bcftoosl_isec')

    bgzipped_tabixed_vcf_file_path_1 = tabix(bgzip(vcf_file_path_1))
    bgzipped_tabixed_vcf_file_path_2 = tabix(bgzip(vcf_file_path_2))

    command = 'bcftools isec {} {} -p {} --threads {}'.format(
        bgzipped_tabixed_vcf_file_path_1, bgzipped_tabixed_vcf_file_path_2,
        output_directory_path, n_jobs)

    run_command(command)

    return output_directory_path


def bcftools_concat(vcf_file_paths, n_jobs):
    """
    Concatenate .vcf files.
    Arguments
        vcf_file_paths (iterable):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_paths[0] + '.bcftools_concat.vcf'

    bgzipped_tabixed_vcf_file_paths = [
        tabix(bgzip(fn)) for fn in vcf_file_paths
    ]

    command = 'bcftools concat -a {} --threads {} > {}'.format(
        ' '.join(bgzipped_tabixed_vcf_file_paths), n_jobs,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def bcftools_rename_chr(vcf_file_path, chromosome_map_file_path, n_jobs=1):
    """
    Rename chromosomes in .vcf file.
    Arguments:
        vcf_file_path (str):
        chromosome_map_file_path (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.bcftools_rename_chr.vcf'

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = 'bcftools annotate --rename-chrs {} --threads {} {} > {}'.format(
        chromosome_map_file_path, n_jobs, bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def bcftools_extract_chromosomes(vcf_file_path,
                                 chromosomes=CHROMOSOMES,
                                 n_jobs=1):
    """
    Extract chromosomes (chromosome 1 to 22, X, Y, and MT) from .vcf file.
    Arguments:
        vcf_file_path (str):
        chromosomes (iterable):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.bcftools_extract_chromosomes.vcf'

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = 'bcftools view -r {} --threads {} {} > {}'.format(
        ','.join(chromosomes), n_jobs, bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def bcftools_filter(vcf_file_path, qual=60, dp=30, n_jobs=1):
    """
    Annotate .vcf file with annotaiton.
    Arguments:
        vcf_file_path (str):
        qual (int):
        dp (int):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.bcftools_filter_qual{}_dp{}.vcf'.format(
        qual, dp)

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = 'bcftoosl view -i \'{}<QUAL & {}<DP\' --threads {} {} > {}'.format(
        qual, dp, n_jobs, bgzipped_tabixed_vcf_file_path, output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def bcftools_annotate(vcf_file_path,
                      annotation_file_path,
                      annotation_arguments='-c ID',
                      n_jobs=1):
    """
    Annotate .vcf file with annotaiton file path.
    Arguments:
        vcf_file_path (str):
        annotation_file_path (str):
        annotation_arguments (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.bcftools_annotate_{}.vcf'.format(
        split(annotation_file_path)[-1])

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = 'bcftools annotate -a {} {} --threads {} {} > {}'.format(
        annotation_file_path, annotation_arguments, n_jobs,
        bgzipped_tabixed_vcf_file_path, output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def snpeff(vcf_file_path, genomic_assembly='GRCh38.86'):
    """
    Annotate .vcf file.
    Arguments:
        vcf_file_path (str):
        genomic_assembly (str): 'GRCh38.86' | 'GRCh37.75'
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.snpeff.vcf'

    if vcf_file_path.endswith('.gz'):
        bgzipped_tabixed_vcf_file_path = tabix(vcf_file_path)
    else:
        bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = 'snpEff -v -noLog -s {}.html {} {} > {}'.format(
        output_vcf_file_path, genomic_assembly, bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def snpsift(vcf_file_path, annotation_file_path):
    """
    Annotate .vcf file with annotaiton file path.
    Arguments:
        vcf_file_path (str):
        annotation_file_path (str):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.snpsift_{}.vcf'.format(
        split(annotation_file_path)[-1])

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = 'SnpSift annotate -noDownload -v -noLog {} {} > {}'.format(
        annotation_file_path, bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))
