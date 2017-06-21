from .helper.helper.subprocess_ import run_command


def concat_snp_indel(snp_file_path, indel_file_path, output_file_path):
    """
    Concatenate SNP.VCF and InDel.VCF.
    :param snp_file_path: str
    :param indel_file_path: str
    :param output_file_path: str
    :return: str
    """

    return concat_sort([snp_file_path, indel_file_path], output_file_path)


def concat_sort(file_paths, output_file_path):
    """
    Concatenate VCFs.
    :param file_paths: iterable
    :param output_file_path: str
    :return: str
    """

    output_file_path = mark_file_path(output_file_path, 'concat_sort',
                                      '.variant')
    file_paths = [bgzip_tabix(fn) for fn in file_paths]

    command = 'bcftools concat -a ' + ' '.join(file_paths) + ' > {}'.format(
        output_file_path)

    run_command(command)

    return bgzip_tabix(output_file_path)


def isec(file_path1, file_path2, output_directory):
    """

    :param file_path1:
    :param file_path2:
    :param output_directory:
    :return: None
    """

    file_path1 = bgzip_tabix(file_path1)
    file_path2 = bgzip_tabix(file_path2)

    command = 'bcftools isec -O z -p {} {} {}'.format(output_directory,
                                                      file_path1, file_path2)

    run_command(command)


def remap_38(file_path, source_assembly):
    """
    Re-map genomic coordinates of <file_path> based on GRCh38 using picard LiftoverVcf.
    :param file_path: str;
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

    output_file_path = mark_file_path(file_path, mark, '.variant')
    reject_file_path = mark_file_path(file_path, '{}_rejected'.format(mark),
                                      '.variant')

    command = '{} LiftoverVcf INPUT={} OUTPUT={} REJECT={} CHAIN={} REFERENCE_SEQUENCE={}'.format(
        PICARD, file_path, output_file_path, reject_file_path, chain,
        path_target_assembly)

    run_command(command)

    bgzip_tabix(reject_file_path)
    return bgzip_tabix(output_file_path)


def rename_chr_sort(file_path):
    """
    Rename chromosomes.
    :param file_path:
    :return: str;
    """

    output_file_path = mark_file_path(file_path, 'rename_chr_sort', '.variant')
    file_path = bgzip_tabix(file_path)

    command = 'bcftools annotate --rename-chrs {} {} -o {}'.format(
        PATH_CHROMOSOME_MAP, file_path, output_file_path)

    run_command(command)

    return bgzip_tabix(output_file_path)


def extract_chr(file_path, chromosome_format):
    """
    Extract chromosome 1 to 22, X, Y, and MT from <file_path>, and bgzip and tabix the extracted VCF.
    :param file_path:
    :param chromosome_format:
    :return: str;
    """

    output_file_path = mark_file_path(file_path, 'extract_chr', '.variant')
    file_path = bgzip_tabix(file_path)

    command_template = 'bcftools view -r {} {} -o {}'
    if chromosome_format == 'chr#':
        command = command_template.format(','.join(CHROMOSOMES_CHR), file_path,
                                          output_file_path)
    elif chromosome_format == '#':
        command = command_template.format(','.join(CHROMOSOMES), file_path,
                                          output_file_path)
    else:
        raise ValueError('Chromosome format {} not found in (chr#, #)'.format(
            chromosome_format))

    run_command(command)

    return bgzip_tabix(output_file_path)


def snpeff(file_path, genomic_assembly):
    """
    Annotate VCF <file_path> using SNPEff.
    :param file_path:
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

    output_file_path = mark_file_path(file_path, 'snpeff', '.variant')

    file_path = bgzip_tabix(file_path)

    command = '{} -noDownload -v -noLog -s {}.html {} {} > {}'.format(
        SNPEFF, output_file_path[:-len('.variant')], genomic_assembly,
        file_path, output_file_path)

    run_command(command)

    return bgzip_tabix(output_file_path)


def snpsift(file_path, annotation):
    """
    Annotate VCF <file_path> using SNPSift.
    :param file_path:
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

    output_file_path = mark_file_path(file_path, mark, '.variant')

    file_path = bgzip_tabix(file_path)

    command = '{} annotate -noDownload -v -noLog {} {} {} > {}'.format(
        SNPSIFT, flag, path_annotation, file_path, output_file_path)

    run_command(command)

    return bgzip_tabix(output_file_path)


def annotate(file_path, genomic_assembly, pipeline):
    """
    Annotate <file_path> VCF with SNPEff and ClinVar (SNPSift).
    :param file_path: str;
    :param genomic_assembly: str;
    :param pipeline: str;
    :return: str;
    """

    if pipeline == 'snpeff':
        print(
            '\n**************************************************************************'
        )
        output_file_path = snpeff(file_path, genomic_assembly)

    elif pipeline == 'snpeff-clinvar':
        print(
            '\n**************************************************************************'
        )
        output_file_path = snpeff(file_path, genomic_assembly)
        print(
            '\n**************************************************************************'
        )
        output_file_path = snpsift(output_file_path, 'clinvar')

    elif pipeline == 'dbsnp-snpeff-clinvar':
        print(
            '\n**************************************************************************'
        )
        output_file_path = snpsift(file_path, 'dbsnp')
        print(
            '\n**************************************************************************'
        )
        output_file_path = snpeff(output_file_path, genomic_assembly)
        print(
            '\n**************************************************************************'
        )
        output_file_path = snpsift(output_file_path, 'clinvar')

    else:
        raise ValueError(
            'Unknown pipeline {}; choose from (snpeff, snpeff-clinvar, dbsnp-snpeff-clinvar).'.
            format(pipeline))

    return bgzip_tabix(output_file_path)
