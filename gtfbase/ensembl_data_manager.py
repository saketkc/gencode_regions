"""
This module will contain methods to fetch data from Ensembl. It will use the FTP module to connect with Ensembl FTP. It will provide
methods to get list of available GTF files and also to download a file from Ensembl FTP.
"""
from ftplib import FTP
from . import utils, consts
import os

class EnsemblDataManager(object):
    """
    """
    def __init__(self, release="current"):
        """
        Constructor for EnsemblDataManager.
        """
        self._release = release
        self._ftp = FTP("ftp.ensembl.org")
        self._ftp.login()

    def go_to_release_directory(self):
        """
        Changes ftp directory to release directory.
        """
        if self._release != "current":
            self._ftp.cwd("/pub/realease-%s/gtf/" % self._release)
        else:
            self._ftp.cwd("/pub/current_gtf/")

    def get_list(self):
        """
        Returns the list of species for the given build from ensembl.
        """
        self.go_to_release_directory()
        return self._ftp.nlst()

    def get_gtf_list(self, species_name):
        """
        Provides the list of different types of gtfs, the species have.
        """
        self.go_to_release_directory()
        self._ftp.cwd(species_name)
        file_list = self._ftp.nlst()
        gtf_files = [file_name for file_name in file_list if
                     file_name.endswith(".gtf.gz")]
        return gtf_files

    def download_gtf_by_filename(self, species_name, file_name, output_dir=None):
        """

        """
        self.go_to_release_directory()
        self._ftp.cwd("%s/" % species_name)
        if output_dir is None:
            output_dir = consts.TEMP_DIR_NAME
        utils.mkdir_p(output_dir)
        file_path = os.path.join(output_dir, file_name)
        with open(file_path, 'wb') as f:
            # Todo: handle case when file name is invalid.
            self._ftp.retrbinary('RETR ' + file_name, f.write)
        return utils.extract_gtf_file(file_path)

    def download_gtf(self, species_name, gtf_type="", output_dir=None):
        """
        Download the desired file with refrence to species name and type.
        If provided file name, it directly downloads it.
        """
        self.go_to_release_directory()
        self._ftp.cwd("%s/" % species_name)
        file_list = self._ftp.nlst()
        gtf_files = [file_name for file_name in file_list if
                     file_name.endswith("%s.gtf.gz" % gtf_type)]
        gtf_files.sort(key=lambda x: len(x))
        # Todo: Case: gtf_files list is empty. Throw error
        our_file = gtf_files[0]
        print(our_file)
        return self.download_gtf_by_filename(species_name, our_file, output_dir)


if __name__ == '__main__':
    our_object = EnsemblDataManager()
    our_object.download_gtf("felis_catus")
