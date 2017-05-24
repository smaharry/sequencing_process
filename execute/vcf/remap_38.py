import argparse

from ...variant import remap_38

parser = argparse.ArgumentParser()
parser.add_argument('input_fname', help='Input compression')
parser.add_argument('source_assembly', help='Source genomic assembly {grch37, hg19}.')
args = parser.parse_args()

remap_38(args.input_fname, args.source_assembly)
