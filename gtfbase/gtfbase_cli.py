import click
from simple_term_menu import TerminalMenu
# from ensembl_data_manager import EnsemblDataManager
# from bed_tool_builder import BedToolBuilderFactory
# import os
# from gene_db import GeneDB
# from consts import TEMP_DIR_NAME
from gtfbase.ensembl_data_manager import EnsemblDataManager
from gtfbase.bed_tool_builder import BedToolBuilderFactory
import os
from gtfbase.gene_db import GeneDB
from gtfbase.consts import TEMP_DIR_NAME



@click.command()
@click.option('--interactive', is_flag=True, help="Please select this if you want"
                                                   "an interactive environment")
@click.option("--path", "-p", type=str, help="path to the gtf file")
@click.option('--build', type=(str, int), help="provide the name of the"
                                               " species and the release from ensembl.\n"
                                               "If you want the current release, put release as 0.\n"
                                               "Either path or build should be passed")
@click.option("--output", "-o", type=str, help="if you want to keep the bed files in any other directory, else not required")
@click.option('--features', '-f',
              type=click.Choice(['exons',
                                 'introns',
                                 'gene',
                                 'CDS',
                                 '5utr',
                                 '3utr',
                                 'start_codon',
                                 'stop_codon',
                                 '3utr_exon',
                                 "3utr_intron",
                                 '5utr_exon',
                                 '5utr_intron',
                                 'all',
                                 'none']), multiple='true')
@click.option('--download_gtf','-dgtf', is_flag=False, help="If you want gtf file to be downloaded, you need to specify build as well.")
def main(interactive, path=None, build=None, output=".", features=None, download_gtf=None):
    """

    """
    if interactive:
        interact()
        return
    click.echo("Please wait your request is processing.")
    if path is not None:
        working(features, path, output)
        return
    name, release = build
    if release is not 0:
        ensembl_object = EnsemblDataManager(release)
    else:
        ensembl_object = EnsemblDataManager()
    gtf_filename = ensembl_object.download_gtf(name, output)
    if download_gtf is True:
        return
    working(features, gtf_filename, output)


def interact():
    release = input("If you do not want latest release please specify the build number else press 'enter':  ")
    if not release.strip():
        ensembl_object = EnsemblDataManager()
    else:
        ensembl_object = EnsemblDataManager(release)
    print("Fetching available species...")
    species_list = ensembl_object.get_species_list()
    species_menu = TerminalMenu(species_list, search_key=None)
    menu_entry_index = species_menu.show()
    print("You have selected %s..." % species_list[menu_entry_index])
    species_name = species_list[menu_entry_index]
    gtf_files = ensembl_object.get_gtf_list(species_name)

    gtf_file_name = gtf_files[0]
    if len(gtf_files) > 1:
        print("Choose one GTF file:")
        gtf_menu = TerminalMenu(gtf_files, search_key=None)
        gtf_entry_idx = gtf_menu.show()
        gtf_file_name = gtf_files[gtf_entry_idx]
    print('%s file will be used...' % gtf_file_name)
    gtf_file_name = ensembl_object.download_gtf_by_filename(species_name, gtf_file_name)
    print("Choose features:")
    feature_list = ['exon',
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

    feature_menu = TerminalMenu(
        feature_list, multi_select=True,
        show_multi_select_hint=True,
    )
    menu_entry_indices = feature_menu.show()
    required_features = [feature_list[x] for x in menu_entry_indices]
    print("You have selected following features: {0}"
          .format(required_features))
    working(required_features, gtf_file_name, ".")


def working(features, gtf, prefix=TEMP_DIR_NAME):
    gene_db = GeneDB(gtf)
    for feature in features:
        # print("Working on {0}".format(feature))
        f_bedtool = BedToolBuilderFactory.get_builder(feature).generate_bedtool(gene_db)
        f_bedtool.saveas(os.path.join(prefix, '%s.bed.gz' % feature))
        # print("Done {0}".format(feature))


if __name__ == '__main__':
    main()



