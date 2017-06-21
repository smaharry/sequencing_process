from .helper.helper.subprocess_ import run_command


def align(hisat2_index_file_path_prefix,
          unpaired_file_path=None,
          paired_forward_file_path=None,
          paired_reverse_file_path=None,
          n_jobs=1):
    """
    Align unpaired .FASTQ or paired .FASTQs.
    :param hisat2_index_file_path_prefix: str
    :param unpaired_file_path: str
    :param paired_forward_file_path: str
    :param paired_reverse_file_path: str
    :param n_jobs: int
    :return: str
    """

    if unpaired_file_path:
        sample_command = '-U {}'.format(unpaired_file_path)
        output_file_path = unpaired_file_path + '.hisat2.sam'

    elif paired_forward_file_path and paired_reverse_file_path:
        sample_command = '-1 {} -2 {}'.format(paired_forward_file_path,
                                              paired_reverse_file_path)
        output_file_path = paired_forward_file_path + '.hisat2.sam'

    else:
        raise ValueError(
            'Need either unpaired sample or paired forward & reverse samples.')

    command = 'hisat2 --dta-cufflinks -p {} -x {} {} -S {}'.format(
        n_jobs, hisat2_index_file_path_prefix, sample_command,
        output_file_path)

    run_command(command)

    return output_file_path
