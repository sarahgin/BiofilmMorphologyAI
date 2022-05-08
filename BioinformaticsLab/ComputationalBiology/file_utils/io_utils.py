import os


def create_dir_if_not_exists(filepath):
    directory = os.path.dirname(filepath)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)