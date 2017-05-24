import argparse

from ...compression import bgzip_tabix

parser = argparse.ArgumentParser()
parser.add_argument('input_fname', help='Input compression')
args = parser.parse_args()

bgzip_tabix(args.input_fname)
