from bed_tool_builder import BedToolBuilderFactory
from gene_db import GeneDB
import os


def working(features, gtf, prefix='.'):
    gene_db = GeneDB(gtf)
    print('gene_db obj done')
    for feature in features:
        f_bedtool = BedToolBuilderFactory.get_builder(feature)._generate_bedtool_helper(gene_db)
        f_bedtool.saveas(os.path.join(prefix, '%s.bed.gz' % feature))


def lets(gtf):
    gene_db = GeneDB(gtf)
    return BedToolBuilderFactory.get_builder('exon')._generate_bedtool_helper(gene_db)


if __name__ == '__main__':
    feature_list = ['exon',
                    'intron',
                    'gene',
                    'CDS',
                    '5utr',
                    '3utr',
                    'start_codon',
                    'stop_codon']
    # working(feature_list, 'Drosophila_melanogaster.BDGP6.32.104.gtf')
    f_bedtool = lets('Drosophila_melanogaster.BDGP6.32.104.gtf')

