"""
This module will contain methods to fetch data from Ensembl. It will use the FTP module to connect with Ensembl FTP. It will provide
methods to get list of available GTF files and also to download a file from Ensembl FTP.
"""
from ftplib import FTP
from gtfbase.utils import mkdir_p, extract_gtf_file, SpeciesNameError
from gtfbase.consts import TEMP_DIR_NAME
import os


class EnsemblDataManager(object):
    """
    This class contains methods to fetch data from Ensembl. It uses the FTP module to connect with Ensembl FTP. It
    provides method to get list of available GTF files and also to download a file from Ensembl FTP.
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

    def get_species_list(self):
        """
        Provides the list of species for the given build from ensembl.

        :return: The list of species
        :rtype: list
        """
        self.go_to_release_directory()
        species_list = self._ftp.nlst()
        return species_list

    def get_gtf_list(self, species_name):
        """
        Provides the list of different types of gtfs currently available with the given species.

        :param species_name: Name of the species of interest
        :type species_name: str
        :return: The list of gtf available
        :rtype: list
        """
        species_name = species_name.lower()
        self.go_to_release_directory()
        try:
            self._ftp.cwd("%s/" % species_name)
        except:
            raise SpeciesNameError("The species name %s is invalid." % species_name)
        file_list = self._ftp.nlst()
        gtf_list = [file_name for file_name in file_list if
                    file_name.endswith(".gtf.gz")]
        return gtf_list

    def download_gtf_by_filename(self, species_name, file_name, output_dir=None):
        """
        Downloads the given file, at the given location. And provides the extracted file name.

        :param species_name: Species Name
        :type species_name: str
        :param file_name: Name of the gtf file
        :type file_name: str
        :param output_dir: Path to the directory, where you want to save the file
        :type output_dir: str
        :return: extracted_file_name: Name of the extracted gtf file
        :rtype: str
        """
        species_name = species_name.lower()
        self.go_to_release_directory()
        try:
            self._ftp.cwd("%s/" % species_name)
        except:
            raise SpeciesNameError("The species name %s is invalid." % species_name)
        if output_dir is None:
            output_dir = TEMP_DIR_NAME
        mkdir_p(output_dir)
        file_path = os.path.join(output_dir, file_name)
        with open(file_path, 'wb') as f:
            self._ftp.retrbinary('RETR ' + file_name, f.write)
        extracted_file_name = extract_gtf_file(file_path)
        return extracted_file_name

    def download_gtf(self, species_name, gtf_type="", output_dir=None):
        """
        Downloads at the given location and provides the name of desired file with reference to species name and type of
        gtf. If you already have the name of the file use method "download_gtf_by_filename()".

        :param species_name: Species Name
        :type species_name: str
        :param gtf_type: The type of the gtf, if you want default don't pass anything. Example if it is "chr" then pass
                         chr
        :type gtf_type: str, optional
        :param output_dir: Path to the directory, where you want to save the file
        :type output_dir: str
        :return: extracted_file_name: Name of the extracted gtf file
        :rtype: str
        """
        species_name = species_name.lower()
        self.go_to_release_directory()
        try:
            self._ftp.cwd("%s/" % species_name)
        except:
            raise SpeciesNameError("The species name %s is invalid." % species_name)
        file_list = self._ftp.nlst()
        gtf_files = [file_name for file_name in file_list if
                     file_name.endswith("%s.gtf.gz" % gtf_type)]
        gtf_files.sort(key=lambda x: len(x))
        our_file = gtf_files[0]
        print(our_file)
        extracted_file_name = self.download_gtf_by_filename(species_name, our_file, output_dir)
        return extracted_file_name


if __name__ == '__main__':
    our_object = EnsemblDataManager()

    help(our_object.get_species_list)
