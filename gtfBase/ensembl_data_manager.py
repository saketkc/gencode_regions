"""
This module will contain methods to fetch data from Ensembl. It will use the FTP module to connect with Ensembl FTP. It will provide
methods to get list of available GTF files and also to download a file from Ensembl FTP.
"""
from ftplib import FTP


class EnsemblDataManager(object):
    """
    """
    def __init__(self, release="current"):
        """

        """
        self._release = release
        self._ftp = FTP("ftp.ensembl.org")
        self._ftp.login()

    def go_to_release_directory(self):
        """

        """
        if self._release != "current":
            self._ftp.cwd("/pub/realease-%s/gtf/" % self._release)
        else:
            self._ftp.cwd("/pub/current_gtf/")

    def get_list(self):
        """

        """
        self.go_to_release_directory()
        return self._ftp.nlst()

    def get_gtf_list(self, species_name):
        """

        """
        self.go_to_release_directory()
        self._ftp.cwd(species_name)
        file_list = self._ftp.nlst()
        gtf_files = [file_name for file_name in file_list if
                     file_name.endswith(".gtf.gz")]
        return gtf_files

    def download_gtf(self, species_name, file_name=None):
        """

        """
        #Todo: type dalo abinitio wali file_name.endswith("%s.gtf.gz" % type)] yeh kar dena
        our_file = file_name
        if our_file is None:
            self.go_to_release_directory()
            self._ftp.cwd("%s/" % species_name)
            file_list = self._ftp.nlst()
            gtf_files = [file_name for file_name in file_list if
                         file_name.endswith(".chr.gtf.gz")]
            our_file = gtf_files[0]
        print(our_file)
        with open(our_file, 'wb') as f:
            self._ftp.retrbinary('RETR ' + our_file, f.write)
        #Todo: think about utils.extract_gtf_file(our_file, "./")


if __name__ == '__main__':
    our_object = EnsemblDataManager()
    our_object.download_gtf("felis_catus")
