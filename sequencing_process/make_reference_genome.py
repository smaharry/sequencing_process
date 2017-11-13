full_name = 'GCA_000001405.15_GRCh38_full_analysis_set.fna.gz'
extra_name = 'hs38DH-extra.fa'
full_plus_extra_name = 'GCA_000001405.15_GRCh38_full_plus_hs38DH-extra_analysis_set.fa'

download(''.format(full_name))
download(''.format(extra_name))
'gzip -dc {0} > {2}; cat {1} >> {2}; gzip {2}'.format(full_name, extra_name,
                                                      full_plus_extra_name)
