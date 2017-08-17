from .gz import bgzip, tabix
from .helper.helper.subprocess_ import run_command

CHROMOSOMES_N = []
CHROMOSOMES_CHRN = []
CHAIN_GRCH37_TO_GRCH38_FILE_PATH = ''
CHAIN_HG19_TO_HG38_FILE_PATH = ''
GRCH38_FILE_PATH = ''
HG38_FILE_PATH = ''
EMSEMBL_VCF_FILE_PATH = ''
CLINVAR_FILE_PATH = ''
CHROMOSOME_MAP_FILE_PATH = ''
PICARD = ''
SNPEFF = ''
SNPSIFT = ''


def concat_sort(vcf_file_paths, output_vcf_file_path):
    """
    Concatenate .vcf files.
    Arguments
        vcf_file_paths (iterable):
        output_vcf_file_path (str):
    Returns:
        str:
    """

    bgzipped_tabixed_vcf_file_paths = [
        tabix(bgzip(fn)) for fn in vcf_file_paths
    ]

    command = 'bcftools concat -a {} -O z -o {}'.format(
        ' '.join(bgzipped_tabixed_vcf_file_paths), output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def isec(vcf_file_path_1, vcf_file_path_2, output_directory_path):
    """
    Get the intersection between 2 .vcf files.
    Arguments:
        vcf_file_path_1 (str):
        vcf_file_path_2 (str):
        output_directory_path (str):
    Returns:
        str:
    """

    bgzipped_tabixed_vcf_file_path_1 = tabix(bgzip(vcf_file_path_1))
    bgzipped_tabixed_vcf_file_path_2 = tabix(bgzip(vcf_file_path_2))

    command = 'bcftools isec {} {} -O z -p {}'.format(
        bgzipped_tabixed_vcf_file_path_1, bgzipped_tabixed_vcf_file_path_2,
        output_directory_path)

    run_command(command)

    return output_directory_path


def remap_38(vcf_file_path, source_assembly):
    """
    Re-map .vcf file coordinates to GRCh38 using picard LiftoverVcf.
    Arguments:
        vcf_file_path (str):
        source_assembly (str): 'hg19' | 'GRCh37'
    Returns:
        str:
    """

    if source_assembly == 'hg19':
        chain = CHAIN_HG19_TO_HG38_FILE_PATH
        target_assembly_file_path = HG38_FILE_PATH

    elif source_assembly == 'GRCh37':
        chain = CHAIN_GRCH37_TO_GRCH38_FILE_PATH
        target_assembly_file_path = GRCH38_FILE_PATH

    else:
        raise ValueError('Unknown source_assembly {}.'.format(source_assembly))

    output_vcf_file_path = vcf_file_path + '.remap_38.vcf'
    rejected_vcf_file_path = vcf_file_path + '.remap_38.rejected.vcf'

    command = '{} LiftoverVcf INPUT={} OUTPUT={} REJECT={} CHAIN={} REFERENCE_SEQUENCE={}'.format(
        PICARD, vcf_file_path, output_vcf_file_path, rejected_vcf_file_path,
        chain, target_assembly_file_path)

    run_command(command)

    tabix(bgzip(rejected_vcf_file_path))

    return tabix(bgzip(output_vcf_file_path))


def rename_chr_sort(vcf_file_path):
    """
    Rename chromosomes in .vcf file.
    Arguments:
        vcf_file_path (str):
    Returns:
        str:
    """

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    output_vcf_file_path = vcf_file_path + '.rename_chr_sort.vcf'

    command = 'bcftools annotate --rename-chrs {} {} -O z -o {}'.format(
        CHROMOSOME_MAP_FILE_PATH, bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def extract_chromosomes_1_to_22_X_Y_MT(vcf_file_path, chromosome_format):
    """
    Extract chromosome 1 to 22, X, Y, and MT from .vcf file.
    Arguments:
        vcf_file_path (str):
        chromosome_format (str): 'chrN' | 'N'
    Returns:
        str:
    """

    output_vcf_file_path = vcf_file_path + '.extract_chromosomes_1_to_22_X_Y_MT.vcf'

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


def snpeff(vcf_file_path, genomic_assembly='GRCh38'):
    """
    Annotate .vcf file using snpeff.
    Arguments:
        vcf_file_path (str):
        genomic_assembly (str): 'GRCh38' | 'GRCh37'
    Returns:
        str:
    """

    if genomic_assembly == 'GRCh37':
        genomic_assembly = 'GRCh37.75'
    elif genomic_assembly == 'GRCh38':
        genomic_assembly = 'GRCh38.82'
    else:
        raise ValueError(
            'Invalid genomic_assembly {}.'.format(genomic_assembly))

    output_vcf_file_path = vcf_file_path + '.snpeff.vcf'

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = '{} -noDownload -v -noLog -s {}.html {} {} > {}'.format(
        SNPEFF, output_vcf_file_path, genomic_assembly,
        bgzipped_tabixed_vcf_file_path, output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))


def snpsift(vcf_file_path, annotation):
    """
    Annotate .vcf file with annotaiton using snpsift.
    Arguments:
        vcf_file_path (str):
        annotation (str): 'ENSEMBL' | 'ClinVar'
    Returns:
        str:
    """

    if annotation == 'ensembl':
        mark = 'ensembl'
        annotation_file_path = EMSEMBL_VCF_FILE_PATH

    elif annotation == 'clinvar':
        mark = 'clinvar'
        annotation_file_path = CLINVAR_FILE_PATH

    else:
        raise ValueError('annotation has to be one of {dbsnp, clinvar}.')

    output_vcf_file_path = vcf_file_path + '.snpsift_{}.vcf'.format(mark)

    bgzipped_tabixed_vcf_file_path = tabix(bgzip(vcf_file_path))

    command = '{} annotate -noDownload -v -noLog {} {} > {}'.format(
        SNPSIFT, annotation_file_path, bgzipped_tabixed_vcf_file_path,
        output_vcf_file_path)

    run_command(command)

    return tabix(bgzip(output_vcf_file_path))
