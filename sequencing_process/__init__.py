from os.path import dirname, join

from .support.support.subprocess_ import run_command

RESOURCE_DIRECTORY_PATH = join(dirname(dirname(__file__)), 'resource')


def print_and_run_command(*args, **kwargs):

    print()
    return run_command(*args, print_command=True, **kwargs)
