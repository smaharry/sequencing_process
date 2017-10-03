from .support.support.subprocess_ import run_command


def align_dna_fastq_with_hisat2(hisat2_index_file_path_prefix,
                                unpaired_fastq_file_path=None,
                                paired_forward_fastq_file_path=None,
                                paired_reverse_fastq_file_path=None,
                                n_jobs=1):
    """
    Align unpaired or paired forward and reverse DNA .fastq files.
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
        output_sam_file_path = unpaired_fastq_file_path + '.hisat2.sam'

    elif paired_forward_fastq_file_path and paired_reverse_fastq_file_path:
        sample_command = '-1 {} -2 {}'.format(paired_forward_fastq_file_path,
                                              paired_reverse_fastq_file_path)
        output_sam_file_path = paired_forward_fastq_file_path + '.hisat2.sam'

    else:
        raise ValueError(
            'Need unpaired or paired forward and reverse DNA .fastq file paths.'
        )

    command = 'hisat2 --no-spliced-alignment {} -x {} -S {} -p {}'.format(
        sample_command, hisat2_index_file_path_prefix, output_sam_file_path,
        n_jobs)

    run_command(command)

    return output_sam_file_path


def align_rna_fastq_with_hisat2(hisat2_index_file_path_prefix,
                                unpaired_fastq_file_path=None,
                                paired_forward_fastq_file_path=None,
                                paired_reverse_fastq_file_path=None,
                                n_jobs=1):
    """
    Align unpaired or paired forward and reverse RNA .fastq files.
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
        output_sam_file_path = unpaired_fastq_file_path + '.hisat2.sam'

    elif paired_forward_fastq_file_path and paired_reverse_fastq_file_path:
        sample_command = '-1 {} -2 {}'.format(paired_forward_fastq_file_path,
                                              paired_reverse_fastq_file_path)
        output_sam_file_path = paired_forward_fastq_file_path + '.hisat2.sam'

    else:
        raise ValueError(
            'Need unpaired or paired forward and reverse RNA .fastq file paths.'
        )

    command = 'hisat2 --dta-cufflinks {} -x {} -S {} -p {}'.format(
        sample_command, hisat2_index_file_path_prefix, output_sam_file_path,
        n_jobs)

    run_command(command)

    return output_sam_file_path
