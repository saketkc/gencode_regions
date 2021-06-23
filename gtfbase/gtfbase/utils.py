"""
This module will contain some general utility functions for our library.
"""
import os
import errno
import tarfile


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def extract_gtf_file(file_name, directory):
    file = tarfile.open(file_name)
    file.extractall(directory)
    file.close()
