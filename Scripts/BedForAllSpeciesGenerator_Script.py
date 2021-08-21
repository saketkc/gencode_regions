from gtfbase import EnsemblDataManager
from gtfbase import GeneDB
from gtfbase import BedToolBuilderFactory
import os
import time
import errno

species_list = [
    "Pan_troglodytes",
    "Macaca_mulatta",
    "Gorilla_gorilla",
    "Pan_paniscus",
    "Callithrix_jacchus",
    "Tupaia_belangeri",
    "Mus_musculus",
    "Rattus_norvegicus",
    "Oryctolagus_cuniculus",
    "Bos_taurus",
    "Sus_scrofa",
    "Gallus_gallus",
    "Xenopus_tropicalis",
    "Anolis_carolinensis",
    "Danio_rerio",
    "Drosophila_melanogaster",
    "Caenorhabditis_elegans"
]
time_taken_by_species = {}

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def script(version):
    path = os.path.join(".", "SpeciesBed") #set a desired path for the directory that will contain all the species
    mkdir_p(path)
    for species in species_list:
      print("Processing species %s" % species)
      path = os.path.join(".", "SpeciesBed", "%s%s" % (species, version))   #
      mkdir_p(path)
      try:
        ensembl_object = EnsemblDataManager()
        gtf = ensembl_object.download_gtf(species)
        features = ['exon',
                    'intron',
                    'gene',
                    'CDS',
                    '5utr',
                    '3utr',
                    'start_codon',
                    'stop_codon',
                    '3utr_exon',
                    "3utr_intron",
                    '5utr_exon',
                    '5utr_intron']
        gene_db = GeneDB(gtf)
        start = time.time()
        for feature in features:
                f_bedtool = BedToolBuilderFactory.get_builder(feature).generate_bedtool(gene_db)
                f_bedtool.saveas(os.path.join(path, '%s.bed.gz' % feature))
        end = time.time()
        time_taken_by_species[species] = (end - start)/60
        print("Time taken by", species, ":", time_taken_by_species[species], "minutes.")
      except Exception as e:
        print(e)
        continue

    for species in species_list:
      print("Time taken by", species, ":", time_taken_by_species[species], "minutes.")


if __name__ == '__main__':
    VERSION = ""
    script(VERSION)
