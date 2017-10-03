from argparse import ArgumentParser

# ==============================================================================
# Set up arguments
# ==============================================================================
parser = ArgumentParser()

parser.add_argument('process', help='Process name.')

parser.add_argument('-vcf', '--vcf-file-path', help='.vcf file path.')

# ==============================================================================
# Parse arguments
# ==============================================================================
args = parser.parse_args()

process = args.process

vcf_file_path = args.vcf_file_path

# ==============================================================================
# Process
# ==============================================================================
