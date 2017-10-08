from . import CHROMOSOMES
from .support.support.subprocess_ import run_command


def concatenate_vcf_gzs_using_bcftools(vcf_gz_file_paths, n_jobs=1):
    """
    Concatenate .vcf.gz files using bcftools.
    Arguments
        vcf_gz_file_paths (iterable):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_gz_file_path = vcf_gz_file_paths[0] + '.concatenate_vcf_gzs_using_bcftools.vcf.gz'

    command = 'bcftools concat -a --threads {} {} | bgzip -cf -@ {} > {}; tabix -f {}'.format(
        n_jobs, ' '.join(vcf_gz_file_paths), n_jobs, output_vcf_gz_file_path,
        output_vcf_gz_file_path)
    run_command(command)

    return output_vcf_gz_file_path


def extract_vcf_gz_chromosomes_using_bcftools(vcf_gz_file_path,
                                              chromosomes=CHROMOSOMES,
                                              n_jobs=1):
    """
    Extract chromosomes from .vcf.gz file using bcftools.
    Arguments:
        vcf_gz_file_path (str):
        chromosomes (iterable):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_gz_file_path = vcf_gz_file_path + '.extract_vcf_gz_chromosomes_using_bcftools.vcf.gz'

    command = 'bcftools view -r {} --threads {} {} | bgzip -cf -@ {} > {}; tabix -f {}'.format(
        ','.join(chromosomes), n_jobs, vcf_gz_file_path, n_jobs,
        output_vcf_gz_file_path, output_vcf_gz_file_path)
    run_command(command)

    return output_vcf_gz_file_path


def filter_vcf_gz_using_bcftools(vcf_gz_file_path, qual=60, dp=30, n_jobs=1):
    """
    Filter .vcf.gz file using bcftools.
    Arguments:
        vcf_gz_file_path (str):
        qual (int):
        dp (int):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_gz_file_path = vcf_gz_file_path + '.filter_vcf_gz_using_bcftools.vcf.gz'

    command = 'bcftoosl view -i \'{}<QUAL & {}<DP\' --threads {} {} | bgzip -cf -@ {} > {}; tabix -f {}'.format(
        qual, dp, n_jobs, vcf_gz_file_path, n_jobs, output_vcf_gz_file_path,
        output_vcf_gz_file_path)
    run_command(command)

    return output_vcf_gz_file_path


def annotate_vcf_gz_using_snpeff(vcf_gz_file_path,
                                 genomic_assembly='GRCh38.90',
                                 maximum_memory='8G',
                                 n_jobs=1):
    """
    Annotate .vcf.gz file using snpeff.
    Arguments:
        vcf_gz_file_path (str):
        genomic_assembly (str):
        maximum_memory (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_gz_file_path = vcf_gz_file_path + '.annotate_vcf_gz_using_snpeff.vcf.gz'

    command = 'snpEff -Xmx{} -s {}.html -v -noLog {} {} | bgzip -cf -@ {} > {}; tabix -f {}'.format(
        maximum_memory, output_vcf_gz_file_path, genomic_assembly,
        vcf_gz_file_path, n_jobs, output_vcf_gz_file_path,
        output_vcf_gz_file_path)
    run_command(command)

    return output_vcf_gz_file_path


def annotate_vcf_gz_using_bcftools(vcf_gz_file_path,
                                   annotation_file_path,
                                   additional_arguments='-c ID',
                                   n_jobs=1):
    """
    Annotate .vcf.gz file using bcftools.
    Arguments:
        vcf_gz_file_path (str):
        annotation_file_path (str):
        additional_arguments (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_vcf_gz_file_path = vcf_gz_file_path + '.annotate_vcf_gz_using_bcftools.vcf.gz'

    command = 'bcftools annotate -a {} --threads {} {} {} | bgzip -cf -@ {} > {}; tabix -f {}'.format(
        annotation_file_path, n_jobs, additional_arguments, vcf_gz_file_path,
        n_jobs, output_vcf_gz_file_path, output_vcf_gz_file_path)
    run_command(command)

    return output_vcf_gz_file_path
