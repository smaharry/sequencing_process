from inspect import stack
from os.path import dirname, isdir, isfile

from ._print_and_run_command import _print_and_run_command
from .bgzip_and_tabix import bgzip_and_tabix
from .process_vcf_gz import concatenate_vcf_gzs_using_bcftools_concat
from .support.support.multiprocess import multiprocess


def sort_and_index_bam_using_samtools_sort_and_index(
        bam_file_path,
        remove_input_bam_file_path=False,
        n_job=1,
        output_bam_file_path=None,
        overwrite=False,
):

    if output_bam_file_path is None:

        output_bam_file_path = '{}/{}.bam'.format(
            dirname(bam_file_path),
            stack()[0][3],
        )

    if not overwrite and isfile(output_bam_file_path):

        raise FileExistsError(output_bam_file_path)

    _print_and_run_command('samtools sort --threads {} {} > {}'.format(
        n_job,
        bam_file_path,
        output_bam_file_path,
    ))

    if remove_input_bam_file_path:

        _print_and_run_command('rm -rf {}'.format(bam_file_path))

    return index_bam_using_samtools_index(
        output_bam_file_path,
        n_job=n_job,
        overwrite=overwrite,
    )


def index_bam_using_samtools_index(
        bam_file_path,
        n_job=1,
        overwrite=False,
):

    output_bai_file_path = bam_file_path + '.bai'

    if not overwrite and isfile(output_bai_file_path):

        raise FileExistsError(output_bai_file_path)

    _print_and_run_command('samtools index -@ {} {}'.format(
        n_job,
        bam_file_path,
    ))

    return bam_file_path


def mark_duplicates_in_bam_using_picard_markduplicates(
        bam_file_path,
        memory='8G',
        remove_duplicates=False,
        remove_input_bam_file_path_and_its_index=False,
        n_job=1,
        output_bam_file_path=None,
        overwrite=False,
):

    if output_bam_file_path is None:

        output_bam_file_path = '{}/{}.bam'.format(
            dirname(bam_file_path),
            stack()[0][3],
        )

    if not overwrite and isfile(output_bam_file_path):

        raise FileExistsError(output_bam_file_path)

    metrics_file_path = output_bam_file_path + '.metrics'

    _print_and_run_command(
        'picard -Xmx{} MarkDuplicates REMOVE_DUPLICATES={} INPUT={} OUTPUT={} METRICS_FILE={}'.
        format(
            memory,
            str(remove_duplicates).lower(),
            bam_file_path,
            output_bam_file_path,
            metrics_file_path,
        ))

    if remove_input_bam_file_path_and_its_index:

        _print_and_run_command('rm -rf {}'.format(bam_file_path))

        _print_and_run_command('rm -rf {}'.format(bam_file_path + '.bai'))

    print('{}:'.format(metrics_file_path))

    with open(metrics_file_path) as metric_file:

        print(metric_file.read())

    return index_bam_using_samtools_index(
        output_bam_file_path,
        n_job=n_job,
        overwrite=overwrite,
    )


def check_bam_using_samtools_flagstat(
        bam_file_path,
        n_job=1,
        output_file_path=None,
        overwrite=False,
):

    if output_file_path is None:

        output_file_path = bam_file_path + '.flagstat'

    if not overwrite and isfile(output_file_path):

        raise FileExistsError(output_file_path)

    _print_and_run_command('samtools flagstat --threads {} {} > {}'.format(
        n_job,
        bam_file_path,
        output_file_path,
    ))

    print('{}:'.format(output_file_path))

    with open(output_file_path) as output_file:

        print(output_file.read())


def check_fastq_gz_or_bam_using_fastqp(
        fastq_gz_or_bam_file_path,
        kmer_length=7,
        overwrite=False,
):

    plot_zip_prefix_path = fastq_gz_or_bam_file_path + '.plot'

    plot_tsv_file_path = plot_zip_prefix_path + '.tsv'

    if not overwrite:

        if isfile(plot_zip_prefix_path + '.zip'):

            raise FileExistsError(plot_zip_prefix_path + '.zip')

        if isfile(plot_tsv_file_path):

            raise FileExistsError(plot_tsv_file_path)

    _print_and_run_command('fastqp --kmer {} --output {} --text {} {}'.format(
        kmer_length,
        plot_zip_prefix_path,
        plot_tsv_file_path,
        fastq_gz_or_bam_file_path,
    ))


def get_variants_from_bam_using_freebayes_and_multiprocess(
        bam_file_path,
        fasta_file_path,
        regions,
        n_job=2,
        output_vcf_file_path=None,
        overwrite=False,
):

    output_vcf_gz_file_path = concatenate_vcf_gzs_using_bcftools_concat(
        multiprocess(
            get_variants_from_bam_using_freebayes,
            ((
                bam_file_path,
                fasta_file_path,
                region,
                1,
                None,
                overwrite,
            ) for region in regions),
            n_job=n_job,
        ),
        remove_input_vcf_gz_file_paths_and_their_indices=True,
        n_job=n_job,
        output_vcf_file_path=output_vcf_file_path,
        overwrite=overwrite,
    )

    return output_vcf_gz_file_path


def get_variants_from_bam_using_freebayes(
        bam_file_path,
        fasta_file_path,
        regions=None,
        n_job=1,
        output_vcf_file_path=None,
        overwrite=False,
):

    if output_vcf_file_path is None:

        output_vcf_file_path = '{}/{}.vcf'.format(
            dirname(bam_file_path),
            stack()[0][3],
        )

    additional_arguments = []

    if regions is not None:

        additional_arguments.append('--region {}'.format(regions))

    if len(additional_arguments):

        output_vcf_file_path = output_vcf_file_path.replace(
            '.vcf',
            '.{}.vcf'.format(' '.join(additional_arguments).replace(
                ' ',
                '_',
            )),
        )

    if not overwrite and isfile(output_vcf_file_path):

        raise FileExistsError(output_vcf_file_path)

    _print_and_run_command('freebayes --fasta-reference {} {} {} > {}'.format(
        fasta_file_path,
        ' '.join(additional_arguments),
        bam_file_path,
        output_vcf_file_path,
    ))

    return bgzip_and_tabix(
        output_vcf_file_path,
        n_job=n_job,
        overwrite=overwrite,
    )


def get_variants_from_bam_using_strelka(
        bam_file_path,
        fasta_file_path,
        output_directory_path,
        n_job=1,
        overwrite=False,
):

    if isdir(output_directory_path):

        if overwrite:

            _print_and_run_command('rm -fr {}'.format(output_directory_path))

        else:

            raise FileExistsError(output_directory_path)

    bash_file_path = '/tmp/strelka.sh'

    with open(
            bash_file_path,
            'w',
    ) as bash_file:

        bash_file.write('source activate sequencing_process_python2.7 &&\n')

        bash_file.write(
            'configureStrelkaGermlineWorkflow.py --bam {} --referenceFasta {} --runDir {} &&\n'.
            format(
                bam_file_path,
                fasta_file_path,
                output_directory_path,
            ))

        bash_file.write('{}/runWorkflow.py --mode local --jobs {}\n'.format(
            output_directory_path,
            n_job,
        ))

    _print_and_run_command('bash {}'.format(bash_file_path))

    stats_file_path = '{}/results/stats/runStats.tsv'.format(
        output_directory_path)

    print('{}:'.format(stats_file_path))

    with open(stats_file_path) as stats_file:

        print(stats_file.read())

    return '{}/results/variants/variants.vcf.gz'.format(output_directory_path)
