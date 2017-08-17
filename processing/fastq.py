from .helper.helper.subprocess_ import run_command


def align(hisat2_index_file_path_prefix,
          unpaired_fastq_file_path=None,
          paired_forward_fastq_file_path=None,
          paired_reverse_fastq_file_path=None,
          n_jobs=1):
    """
    Align unpaired .fastq file or paired forward and reverse .fastq files.
    Arguments:
        hisat2_index_file_path_prefix (str):
        unpaired_fastq_file_path (str):
        paired_forward_fastq_file_path (str):
        paired_reverse_fastq_file_path (str):
        n_jobs (int):
    Returns:
        str:
    """

    if unpaired_fastq_file_path:
        sample_command = '-U {}'.format(unpaired_fastq_file_path)
        output_file_path = unpaired_fastq_file_path + '.hisat2.sam'

    elif paired_forward_fastq_file_path and paired_reverse_fastq_file_path:
        sample_command = '-1 {} -2 {}'.format(paired_forward_fastq_file_path,
                                              paired_reverse_fastq_file_path)
        output_file_path = paired_forward_fastq_file_path + '.hisat2.sam'

    else:
        raise ValueError(
            'Need unpaired .fastq file path or paired forward and reverse .fastq file paths.'
        )

    command = 'hisat2 --dta-cufflinks -p {} -x {} {} -S {}'.format(
        n_jobs, hisat2_index_file_path_prefix, sample_command,
        output_file_path)

    run_command(command)

    return output_file_path
