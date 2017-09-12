from .subprocess_ import run_command


def shutdown():
    """
    Shutdown.
    Arguments:
        None
    Retuns:
        None
    """

    run_command('sudo shutdown -h now')


def restart():
    """
    Restart.
    Arguments:
        None
    Retuns:
        None
    """

    run_command('sudo reboot')
