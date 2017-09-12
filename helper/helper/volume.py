from os.path import isdir

from .subprocess_ import run_command


def get_volume_name(volume_label):
    """
    Get volume name.
    Arguments:
        volume_label (str):
    Retuns:
        str:
    """

    volume_dict = make_volume_dict()

    for volume_name, d in volume_dict.items():
        if 'LABEL' in d and d['LABEL'] == volume_label:
            return volume_name


def make_volume_dict():
    """
    Make volume dict.
    Arguments:
        None
    Retuns:
        dict: {'/dev/sda1' : {'LABEL': 'MyUSB', 'TYPE': 'ntfs', ...}, ...}
    """

    volume_dict = {}

    for line in run_command('sudo blkid').stdout.strip('\n').split('\n'):

        line = line.split()

        volume_name = line[0][:-1]

        volume_dict[volume_name] = {}

        for fv in line[1:]:

            f, v = fv.replace('\"', '').split('=')

            volume_dict[volume_name][f] = v

    return volume_dict


def mount(volume_name, mount_directory_path):
    """
    Mount the volume_name onto mount_directory_path.
    Arguments:
        volume_name (str):
        mount_directory_path (str):
    Retuns:
        None
    """

    if not isdir(mount_directory_path):
        raise ValueError('{} does not exist. Make it by\n$ sudo mkdir -pv {}'.
                         format(mount_directory_path))

    run_command('sudo mount {} {}'.format(volume_name, mount_directory_path))


def unmount(volume_name=None, mount_directory_path=None):
    """
    Unmount the volume_name.
    Arguments:
        volume_name (str):
        mount_directory_path (str):
    Retuns:
        None
    """

    if volume_name:
        run_command('sudo umount {}'.format(volume_name))

    elif mount_directory_path:
        run_command('sudo umount {}'.format(mount_directory_path))

    else:
        raise ValueError('Need either volume_name or mount_directory_path.')
