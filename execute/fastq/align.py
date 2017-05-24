import argparse

from ...sequence import align

parser = argparse.ArgumentParser()
parser.add_argument('sample', nargs='+', help='Filepath to either unpaired sample (1 filepath) or forward & reverse paired samples (2 filepaths)')
parser.add_argument('--n_threads', dest='n_threads', default=1, help='Number of threads')

args = parser.parse_args()

align(filepaths_sample=args.sample, n_threads=args.n_threads)