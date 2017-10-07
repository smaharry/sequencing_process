from .support.support.subprocess_ import run_command


def align_fastq_gz_with_hisat2(hisat2_index_file_path_prefix,
                               sequence_type,
                               unpaired_fastq_gz_file_path=None,
                               paired_forward_fastq_gz_file_path=None,
                               paired_reverse_fastq_gz_file_path=None,
                               n_jobs=1):
    """
    Align unpaired or paired .fastq.gz files.
    Arguments:
        hisat2_index_file_path_prefix (str):
        sequence_type (str): 'DNA' | 'RNA'
        unpaired_fastq_gz_file_path (str):
        paired_forward_fastq_gz_file_path (str):
        paired_reverse_fastq_gz_file_path (str):
        n_jobs (int):
    Returns:
        str:
    """

    extension = '.hisat2.sort.bam'

    if unpaired_fastq_gz_file_path:
        sample_command = '-U {}'.format(unpaired_fastq_gz_file_path)
        output_bam_file_path = unpaired_fastq_gz_file_path + extension

    elif paired_forward_fastq_gz_file_path and paired_reverse_fastq_gz_file_path:
        sample_command = '-1 {} -2 {}'.format(
            paired_forward_fastq_gz_file_path,
            paired_reverse_fastq_gz_file_path)
        output_bam_file_path = paired_forward_fastq_gz_file_path + extension

    else:
        raise ValueError('Need unpaired or paired .fastq.gz file paths.')

    if sequence_type == 'DNA':
        additional_arguments = '--dta-cufflinks'

    elif sequence_type == 'RNA':
        additional_arguments = '--no-spliced-alignment'

    command = 'hisat2 {} -x {} -p {} {} | samtools view -Sb -@ {} | samtools sort -@ {} > {} & samtools index -@ {} {}'.format(
        sample_command,
        hisat2_index_file_path_prefix,
        n_jobs,
        additional_arguments,
        n_jobs,
        n_jobs,
        output_bam_file_path,
        n_jobs,
        output_bam_file_path, )

    run_command(command)

    return output_bam_file_path
