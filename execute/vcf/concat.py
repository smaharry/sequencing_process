import argparse

from ...variant import concat_sort


parser = argparse.ArgumentParser()
parser.add_argument('vcf_1', help='1st Input VCF compression {.variant, .variant.gz}')
parser.add_argument('vcf_2', help='2nd Input VCF compression {.variant, .variant.gz}')
parser.add_argument('output_fname', help='Output compression')
args = parser.parse_args()

concat_sort([args.vcf_1, args.vcf_2], args.output_fname)
