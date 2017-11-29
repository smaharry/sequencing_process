from os.path import dirname, join

from .support.support.path import clean_path

RESOURCE_DIRECTORY_PATH = clean_path(
    join(dirname(dirname(__file__)), 'resource'))
