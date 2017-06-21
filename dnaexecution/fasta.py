from .helper.helper.subprocess_ import run_command


def get_sequence(file_path, chromosome, start, end):
    """
    Return genomic sequences from region chromosome:start-end.
    :param file_path: str; .FASTA or .FA
    :param chromosome: int | str; chromosome
    :param start: int | str; start position
    :param end: int | str; end position; must be greater than or equal to the
        start position
    :return: str; genomic sequences from region chromosome:start-end
    """

    if int(end) <= int(start):
        raise ValueError(
            'End position {} must be greater than start position {}.'.format(
                end, start))

    command = 'samtools faidx {} {}:{}-{}'.format(file_path, chromosome, start,
                                                  end)

    stdout = run_command(command).stdout

    return ''.join(stdout.split('\n')[1:])
