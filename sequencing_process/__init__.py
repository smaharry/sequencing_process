from os.path import dirname, join

from .support.support.subprocess_ import run_command

RESOURCE_DIRECTORY_PATH = join(dirname(dirname(__file__)), 'resource')


def print_and_run_command(*args, **kwargs):
    """Print and run command."""

    print()
    print('=' * 80)
    return run_command(*args, print_command=True, **kwargs)
    print('=' * 80)
    print()
