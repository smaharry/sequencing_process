from inspect import stack
from os.path import dirname, isfile, join

from . import RESOURCE_DIRECTORY_PATH
from .bgzip_and_tabix import bgzip_and_tabix
from .print_and_run_command import print_and_run_command


def concatenate_vcf_gzs_using_bcftools_concat(
        vcf_gz_file_paths,
        remove_input_vcf_gz_file_paths_and_their_indices=False,
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False):

    if output_vcf_file_path is None:

        output_vcf_file_path = join(
            dirname(vcf_gz_file_paths[0]),
            stack()[0][3] + '.vcf')

    if not overwrite and isfile(output_vcf_file_path + '.gz'):

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


def rename_chromosome_of_vcf_gz_using_bcftools_annotate(
        vcf_gz_file_path,
        map_file_path=join(RESOURCE_DIRECTORY_PATH, 'chrn_n.tsv'),
        remove_input_vcf_gz_file_path_and_its_index=False,
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False):

    if output_vcf_file_path is None:

        output_vcf_file_path = join(
            dirname(vcf_gz_file_path),
            stack()[0][3] + '.vcf')

    if not overwrite and isfile(output_vcf_file_path + '.gz'):

        raise FileExistsError(output_vcf_file_path + '.gz')

    print_and_run_command(
        'bcftools annotate --rename-chrs {} --threads {} {} > {}'.format(
            map_file_path, n_job, vcf_gz_file_path, output_vcf_file_path))

    if remove_input_vcf_gz_file_path_and_its_index:

        print_and_run_command('rm -rf {}'.format(vcf_gz_file_path))

        print_and_run_command('rm -rf {}'.format(vcf_gz_file_path + '.tbi'))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)


def annotate_vcf_gz_using_snpeff(
        vcf_gz_file_path,
        genomic_assembly,
        memory='8G',
        remove_input_vcf_gz_file_path_and_its_index=False,
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False):

    if output_vcf_file_path is None:

        output_vcf_file_path = join(
            dirname(vcf_gz_file_path),
            stack()[0][3] + '.vcf')

    if not overwrite and isfile(output_vcf_file_path + '.gz'):

        raise FileExistsError(output_vcf_file_path + '.gz')

    print_and_run_command(
        'snpEff -Xmx{} -htmlStats {}.stats.html -csvStats {}.stats.csv -t -verbose -noLog {} {} > {}'.
        format(memory, output_vcf_file_path, output_vcf_file_path,
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

    if output_vcf_file_path is None:

        output_vcf_file_path = join(
            dirname(vcf_gz_file_path),
            stack()[0][3] + '.vcf')

    if not overwrite and isfile(output_vcf_file_path + '.gz'):

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


def filter_vcf_gz_using_bcftools_view(vcf_gz_file_path,
                                      regions=None,
                                      keep_filters=None,
                                      include_expression=None,
                                      n_job=1,
                                      output_vcf_file_path=None,
                                      overwrite=False):

    if output_vcf_file_path is None:

        output_vcf_file_path = join(
            dirname(vcf_gz_file_path),
            stack()[0][3] + '.vcf')

    if not overwrite and isfile(output_vcf_file_path + '.gz'):

        raise FileExistsError(output_vcf_file_path + '.gz')

    additional_arguments = []

    if regions is not None:

        additional_arguments.append('--regions {}'.format(','.join(regions)))

    if keep_filters is not None:

        additional_arguments.append('--apply-filters {}'.format(
            ','.join(keep_filters)))

    if include_expression is not None:

        additional_arguments.append(
            '--include \'{}\''.format(include_expression))

    print_and_run_command('bcftools view {} --threads {} {} > {}'.format(
        ' '.join(additional_arguments), n_job, vcf_gz_file_path,
        output_vcf_file_path))

    return bgzip_and_tabix(
        output_vcf_file_path, n_job=n_job, overwrite=overwrite)
