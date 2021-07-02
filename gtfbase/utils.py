"""
This module will contain some general utility functions for our library.
"""
import os
import errno
import gzip
import shutil

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def extract_gtf_file(file_name):
    """
    Extracts the .gz file and returns the name of the extracted file.
    """
    with gzip.open(file_name, 'rb') as f_in:
        file_name = file_name[:-3]
        with open(file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return file_name


if __name__ == '__main__':
    print(extract_gtf_file("Felis_catus.Felis_catus_9.0.104.gtf.gz"))
