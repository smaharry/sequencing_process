import gzip
from os import remove, rename
from os.path import isfile

from Bio import bgzf


# ======================================================================================================================
# Compression functions
# ======================================================================================================================
def bgzip_tabix(fname):
    """
    bgzip and tabix <fname>.
    :param fname:
    :return: str;
    """

    if fname.endswith('.gz'):
        if not isfile('{}.tbi'.format(fname)):  # tabixed
            tabix(fname)
        return fname
    else:
        bgzip_output_fname = bgzip(fname)
        tabix(bgzip_output_fname)

    return bgzip_output_fname


def bgzip(fname):
    """
    bgzip <fname>.
    :param fname: str;
    :return: str;
    """

    cmd = 'bgzip -f {}'.format(fname)

    run_command(cmd)

    return '{}.gz'.format(fname)


def tabix(fname):
    """
    tabix <fname>.
    :param fname: str;
    :return: str;
    """

    cmd = 'tabix -f {}'.format(fname)

    run_command(cmd)

    return '{}.tbi'.format(fname)


def convert_gzipped_to_bgzipped(gzipped_filename):
    """

    :param gzipped_filename:
    :return: None
    """

    temp_filename = '{}.convert_gzipped_to_bgzipped'.format(gzipped_filename)
    with gzip.open(gzipped_filename, 'rt') as gzipped_file:
        with bgzf.open(temp_filename, 'wt') as temp_bgzipped_file:
            for line in gzipped_file:
                temp_bgzipped_file.write(line)
    remove(gzipped_filename)
    rename(temp_filename, gzipped_filename)


# ======================================================================================================================
# .FASTA or .FA
# ======================================================================================================================
def get_sequence(filepath, chromosome, start, end):
    """
    Return genomic sequences from region specified by chromosome:start-stop in filepath.
    :param filepath: str;
    :param chromosome: int or str; chromosome
    :param start: int or str; start position
    :param end: int or str; end position; must be greater than or equal to start position
    :return: str; genomic sequences from region specified by chromosome:start-stop in filepath
    """

    if start >= end:
        raise ValueError(
            'Starting genomic position must be greater than the ending genomic position.'
        )

    cmd = 'samtools faidx {} {}:{}-{}'.format(filepath, chromosome, start, end)

    stdout = run_command(cmd).stdout

    s = ''.join(stdout.split('\n')[1:])

    return s


# ======================================================================================================================
# .FASTQ
# ======================================================================================================================
def align(filepaths_sample, n_threads=1):
    """

    :param filepaths_sample: str;
    :param path_hisat2_index: str;
    :param n_threads: int;
    :return: None
    """

    if len(filepaths_sample) == 1:
        sample_command = '-U {}'.format(*filepaths_sample)
    elif len(filepaths_sample) == 2:
        sample_command = '-1 {} -2 {}'.format(*filepaths_sample)
    else:
        raise ValueError(
            'Accept either unpaired sample (1 filepath) or forward & reverse samples (2 filepaths).'
        )

    output_filepath = filepaths_sample[0] + '.hisat2.sam'

    command = 'hisat2 --dta-cufflinks -p {} -x {} {} -S {}'.format(
        n_threads, PATH_HISAT2_INDEX, sample_command, output_filepath)
    run_command(command)


# ======================================================================================================================
# .VCF(.GZ)
# ======================================================================================================================
def concat_snp_indel(snp_filepath, indel_filepath, output_fname):
    """
    Concatenate SNP and InDel using BCFTools concat and sort.
    :param snp_filepath: str;
    :param indel_filepath: str;
    :param output_fname: str;
    :return: str;
    """

    return concat_sort([snp_filepath, indel_filepath], output_fname)


def concat_sort(fnames, output_fname):
    """
    Concatenate VCFs <fnames> using BCFTools concat and sort.
    :param fnames:
    :param output_fname: str;
    :return: str;
    """

    output_fname = mark_filename(output_fname, 'concat_sort', '.variant')
    fnames = [bgzip_tabix(fn) for fn in fnames]

    cmd = 'bcftools concat -a ' + ' '.join(fnames) + ' > {}'.format(
        output_fname)

    run_command(cmd)

    return bgzip_tabix(output_fname)


def isec(fname1, fname2, output_directory):
    """

    :param fname1:
    :param fname2:
    :param output_directory:
    :return: None
    """

    fname1 = bgzip_tabix(fname1)
    fname2 = bgzip_tabix(fname2)

    cmd = 'bcftools isec -O z -p {} {} {}'.format(output_directory, fname1,
                                                  fname2)

    run_command(cmd)


def remap_38(fname, source_assembly):
    """
    Re-map genomic coordinates of <fname> based on GRCh38 using picard LiftoverVcf.
    :param fname: str;
    :param source_assembly: str; {hg19, grch37}
    :return: str;
    """

    if source_assembly == 'hg19':
        chain = PATH_CHAIN_HG19_TO_HG38
        path_target_assembly = PATH_HG38
    elif source_assembly == 'grch37':
        chain = PATH_CHAIN_GRCH37_TO_GRCH38
        path_target_assembly = PATH_GRCH38
    else:
        raise ValueError('Unknown assembly_from {}.'.format(source_assembly))

    mark = '{}grch38'.format(source_assembly)

    output_fname = mark_filename(fname, mark, '.variant')
    reject_fname = mark_filename(fname, '{}_rejected'.format(mark), '.variant')

    cmd = '{} LiftoverVcf INPUT={} OUTPUT={} REJECT={} CHAIN={} REFERENCE_SEQUENCE={}'.format(
        PICARD, fname, output_fname, reject_fname, chain, path_target_assembly)

    run_command(cmd)

    bgzip_tabix(reject_fname)
    return bgzip_tabix(output_fname)


def rename_chr_sort(fname):
    """
    Rename chromosomes.
    :param fname:
    :return: str;
    """

    output_fname = mark_filename(fname, 'rename_chr_sort', '.variant')
    fname = bgzip_tabix(fname)

    cmd = 'bcftools annotate --rename-chrs {} {} -o {}'.format(
        PATH_CHROMOSOME_MAP, fname, output_fname)

    run_command(cmd)

    return bgzip_tabix(output_fname)


def extract_chr(fname, chromosome_format):
    """
    Extract chromosome 1 to 22, X, Y, and MT from <fname>, and bgzip and tabix the extracted VCF.
    :param fname:
    :param chromosome_format:
    :return: str;
    """

    output_fname = mark_filename(fname, 'extract_chr', '.variant')
    fname = bgzip_tabix(fname)

    cmd_template = 'bcftools view -r {} {} -o {}'
    if chromosome_format == 'chr#':
        cmd = cmd_template.format(','.join(CHROMOSOMES_CHR), fname,
                                  output_fname)
    elif chromosome_format == '#':
        cmd = cmd_template.format(','.join(CHROMOSOMES), fname, output_fname)
    else:
        raise ValueError('Chromosome format {} not found in (chr#, #)'.format(
            chromosome_format))

    run_command(cmd)

    return bgzip_tabix(output_fname)


def snpeff(fname, genomic_assembly):
    """
    Annotate VCF <fname> using SNPEff.
    :param fname:
    :param genomic_assembly: str;
    :return: str;
    """

    if genomic_assembly == 'grch37':
        genomic_assembly = 'GRCh37.75'
    elif genomic_assembly == 'grch38':
        genomic_assembly = 'GRCh38.82'
    else:
        raise ValueError(
            'Unknown genomic_assembly {}; choose from (grch37, grch38).'.
            format(genomic_assembly))

    output_fname = mark_filename(fname, 'snpeff', '.variant')

    fname = bgzip_tabix(fname)

    cmd = '{} -noDownload -v -noLog -s {}.html {} {} > {}'.format(
        SNPEFF, output_fname[:-len('.variant')], genomic_assembly, fname,
        output_fname)

    run_command(cmd)

    return bgzip_tabix(output_fname)


def snpsift(fname, annotation):
    """
    Annotate VCF <fname> using SNPSift.
    :param fname:
    :param annotation: str; {dpsnp, clinvar}
    :return: str;
    """

    if annotation == 'dbsnp':
        mark = 'dbsnp'
        path_annotation = PATH_DBSNP
        flag = '-noInfo'

    elif annotation == 'clinvar':
        mark = 'clinvar'
        path_annotation = PATH_CLINVAR
        flag = ''

    else:
        raise ValueError('annotation has to be one of {dbsnp, clinvar}.')

    output_fname = mark_filename(fname, mark, '.variant')

    fname = bgzip_tabix(fname)

    cmd = '{} annotate -noDownload -v -noLog {} {} {} > {}'.format(
        SNPSIFT, flag, path_annotation, fname, output_fname)

    run_command(cmd)

    return bgzip_tabix(output_fname)


def annotate(fname, genomic_assembly, pipeline):
    """
    Annotate <fname> VCF with SNPEff and ClinVar (SNPSift).
    :param fname: str;
    :param genomic_assembly: str;
    :param pipeline: str;
    :return: str;
    """

    if pipeline == 'snpeff':
        print(
            '\n**************************************************************************'
        )
        output_fname = snpeff(fname, genomic_assembly)

    elif pipeline == 'snpeff-clinvar':
        print(
            '\n**************************************************************************'
        )
        output_fname = snpeff(fname, genomic_assembly)
        print(
            '\n**************************************************************************'
        )
        output_fname = snpsift(output_fname, 'clinvar')

    elif pipeline == 'dbsnp-snpeff-clinvar':
        print(
            '\n**************************************************************************'
        )
        output_fname = snpsift(fname, 'dbsnp')
        print(
            '\n**************************************************************************'
        )
        output_fname = snpeff(output_fname, genomic_assembly)
        print(
            '\n**************************************************************************'
        )
        output_fname = snpsift(output_fname, 'clinvar')

    else:
        raise ValueError(
            'Unknown pipeline {}; choose from (snpeff, snpeff-clinvar, dbsnp-snpeff-clinvar).'.
            format(pipeline))

    return bgzip_tabix(output_fname)
