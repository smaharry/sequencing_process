from .support.support.subprocess_ import run_command


def _print_and_run_command(*args, **kwargs):

    print()

    return run_command(*args, print_command=True, **kwargs)
