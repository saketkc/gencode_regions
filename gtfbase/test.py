from bed_tool_builder import BedToolBuilderFactory
from gene_db import GeneDB
import os


def working(features, gtf, prefix='.'):
    gene_db = GeneDB(gtf)
    print('gene_db obj done')
    for feature in features:
        f_bedtool = BedToolBuilderFactory.get_builder(feature).generate_bedtool(gene_db)
        f_bedtool.saveas(os.path.join(prefix, '%s.bed.gz' % feature))


if __name__ == '__main__':
    working(['exon'], 'Felis_catus.Felis_catus_9.0.104.chr.gtf')