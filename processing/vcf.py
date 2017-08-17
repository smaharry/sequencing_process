from os.path import join, split

from . import CHROMOSOMES_CHRN, CHROMOSOMES_N
from .gz import bgzip, tabix
from .helper.helper.subprocess_ import run_command

PICARD = ''
SNPEFF = ''
SNPSIFT = ''


def picard_liftovervcf(vcf_file_path, source_assembly,
                       assembly_chain_file_path, target_assembly_file_path):
    """
    Re-map .vcf file coordinates.
    Arguments:
        vcf_file_path (str):
        source_assembly (str): 'hg19' | 'GRCh37'
        assembly_chain_file_path (str):
        target_assembly_file_path (str):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.picard_liftovervcf.vcf'

    rejected_vcf_file_path = vcf_file_path + '.picard_liftovervcf.rejected.vcf'

    command = '{} LiftoverVcf INPUT={} OUTPUT={} REJECT={} CHAIN={} REFERENCE_SEQUENCE={}'.format(
        PICARD, vcf_file_path, output_vcf_file_path, rejected_vcf_file_path,
        assembly_chain_file_path, target_assembly_file_path)

    run_command(command)

    tabix(bgzip(rejected_vcf_file_path))

    return tabix(bgzip(output_vcf_file_path))


def bcftools_isec(vcf_file_path_1, vcf_file_path_2):
    """
    Get the intersection between 2 .vcf files.
    Arguments:
        vcf_file_path_1 (str):
        vcf_file_path_2 (str):
    Returns:
        str:
    """

    output_directory_path = join(split(vcf_file_path_1)[0], 'bcftoosl_isec')

    bgzipped_tabixed_vcf_file_path_1 = tabix(bgzip(vcf_file_path_1))
    bgzipped_tabixed_vcf_file_path_2 = tabix(bgzip(vcf_file_path_2))

    command = 'bcftools isec {} {} -O z -p {}'.format(
        bgzipped_tabixed_vcf_file_path_1, bgzipped_tabixed_vcf_file_path_2,
        output_directory_path)

    run_command(command)

    return output_directory_path


def bcftools_concat(vcf_file_paths):
    """
    Concatenate .vcf files.
    Arguments
        vcf_file_paths (iterable):
        output_vcf_file_path (str):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_paths[0] + '.bcftools_concat.vcf'

    bgzipped_tabixed_vcf_file_paths = [
        tabix(bgzip(fn)) for fn in vcf_file_paths
    ]

    command = 'bcftools concat -a {} -O z -o {}'.format(
        ' '.join(bgzipped_tabixed_vcf_file_paths), output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def bcftools_rename_chr(vcf_file_path, chromosome_map_file_path):
    """
    Rename chromosomes in .vcf file.
    Arguments:
        vcf_file_path (str):
        chromosome_map_file_path (str):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.bcftools_rename_chr.vcf'

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = 'bcftools annotate --rename-chrs {} {} -O z -o {}'.format(
        chromosome_map_file_path, bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def bcftools_extract_chromosomes(vcf_file_path, chromosome_format):
    """
    Extract chromosomes (chromosome 1 to 22, X, Y, and MT) from .vcf file.
    Arguments:
        vcf_file_path (str):
        chromosome_format (str): 'chrN' | 'N'
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.bcftools_extract_chromosomes.vcf'

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    if chromosome_format == 'chrN':
        chromosomes = CHROMOSOMES_CHRN
    elif chromosome_format == 'N':
        chromosomes = CHROMOSOMES_N
    else:
        raise ValueError(
            'Invalid chromosome_format {}'.format(chromosome_format))

    command = 'bcftools view -r {} {} -O z -o {}'.format(
        ','.join(chromosomes), bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def bcftools_filter(vcf_file_path, qual=60, dp=30):
    """
    Annotate .vcf file with annotaiton.
    Arguments:
        vcf_file_path (str):
        qual (int):
        dp (int):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.bcftools_filter_qual{}_dp{}.vcf'.format(
        qual, dp)

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = 'bcftoosl view '.format(bgzipped_tabixed_vcf_file_path, qual, dp,
                                      output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def snpeff(vcf_file_path, genomic_assembly='GRCh38.82'):
    """
    Annotate .vcf file.
    Arguments:
        vcf_file_path (str):
        genomic_assembly (str): 'GRCh38.82' | 'GRCh37.75'
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.snpeff.vcf'

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = '{} -noDownload -v -noLog -s {}.html {} {} > {}'.format(
        SNPEFF, output_vcf_file_path, genomic_assembly,
        bgzipped_tabixed_vcf_file_path, output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def snpsift(vcf_file_path, annotation_file_path):
    """
    Annotate .vcf file with annotaiton.
    Arguments:
        vcf_file_path (str):
        annotation_file_path (str):
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.snpsift_{}.vcf'.format(
        annotation_file_path)

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = '{} annotate -noDownload -v -noLog {} {} > {}'.format(
        SNPSIFT, annotation_file_path, bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))
