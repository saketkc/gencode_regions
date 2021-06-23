"""
This module will contain methods to generate a BedTool object from a given
Gene Dictionary.
"""
import pybedtools
import re


class AbstractBedToolBuilder(object):
    _feature_name = None
    
    def _make_bed(self, gene_db):
        # Todo: if self._feature_name is None: error dena ha

        f_bed = ""
        for gene_id in gene_db.get_gene_list():
            f_regions = []
            for feature in gene_db.gene_dict[gene_id].keys():
                if feature == 'gene':
                    continue
                f = list(gene_db.gene_dict[gene_id][feature][self._feature_name])
                f_regions += f

            merged_f = self._merge_regions(gene_db.feature_db, f_regions)
            renamed_f = self._rename_regions(merged_f, gene_id)
            f_bed += self._create_bed(renamed_f)

        return f_bed

    def generate_bedtool(self, gene_db):
        f_bed = self._make_bed(gene_db)
        f_bedtool = pybedtools.BedTool(f_bed, from_string=True)
        f_bedtool = f_bedtool.remove_invalid().sort()
        return f_bedtool

    def _create_bed(self, regions, bedtype='0'):
        '''Create bed from list of regions
        bedtype: 0 or 1
            0-Based or 1-based coordinate of the BED
        '''
        bedstr = ''
        for region in regions:
            assert len(region.attributes['gene_id']) == 1
            ## GTF start is 1-based, so shift by one while writing
            ## to 0-based BED format
            if bedtype == '0':
                start = region.start - 1
            else:
                start = region.start
            bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(region.chrom,
                                                        start,
                                                        region.stop,
                                                        re.sub('\.\d+', '', region.attributes['gene_id'][0]),
                                                        '.',
                                                        region.strand)
        return bedstr

    def _rename_regions(self, regions, gene_id):
        regions = list(regions)
        if len(regions) == 0:
            return []
        for region in regions:
            region.attributes['gene_id'] = gene_id
        return regions

    def _merge_regions(self, db, regions):
        if len(regions) == 0:
            return []
        merged = db.merge(sorted(list(regions), key=lambda x: x.start))
        return merged

    def _merge_regions_nostrand(self, db, regions):
        if len(regions) == 0:
            return []
        merged = db.merge(sorted(list(regions), key=lambda x: x.start), ignore_strand=True)
        return merged


class GeneBedToolBuilder(AbstractBedToolBuilder):
    _feature_name = 'gene'

    def _make_bed(self, gene_db):
        gene_list = []
        for gene_id in gene_db.get_gene_list():
            gene_list.append(gene_db.gene_dict[gene_id]['gene'])

        gene_bed = self._create_bed(gene_list)

        return gene_bed


class CdsBedToolBuilder(AbstractBedToolBuilder):
    _feature_name = 'CDS'


class ExonBedToolBuilder(AbstractBedToolBuilder):
    _feature_name = 'exon'


class IntronBedToolBuilder(AbstractBedToolBuilder):
    _feature_name = 'exon'

    def _make_bed(self, gene_db):
        intron_bed = ''

        for gene_id in gene_db.get_gene_list():
            f_regions = []
            intron_regions = []
            for feature in gene_db.gene_dict[gene_id].keys():
                if feature == 'gene':
                    continue
                f = list(gene_db.gene_dict[gene_id][feature][self._feature_name])
                merged_f = self._merge_regions(gene_db.feature_db, f)
                introns = gene_db.feature_db.interfeatures(merged_f)
                f_regions += f
                intron_regions += introns

            merged_introns = self._merge_regions(gene_db.feature_db, intron_regions)
            renamed_introns = self._rename_regions(merged_introns, gene_id)

            intron_bed += self._create_bed(renamed_introns)

        return intron_bed


class FivePrimeUtrBedToolBuilder(AbstractBedToolBuilder):
    _feature_name = '5utr'

    def generate_bedtool(self, gene_db):
        utr5_bed = self._make_bed(gene_db)
        utr5_bedtool = pybedtools.BedTool(utr5_bed, from_string=True)
        cds_bedtool = CdsBedToolBuilder().generate_bedtool(gene_db)
        utr5_bedtool = utr5_bedtool.subtract(cds_bedtool).remove_invalid().sort()
        return utr5_bedtool


class ThreePrimeUtrCdsBedToolBuilder(AbstractBedToolBuilder):
    _feature_name = '3utr'

    def generate_bedtool(self, gene_db):
        utr3_bed = self._make_bed(gene_db)
        utr3_bedtool = pybedtools.BedTool(utr3_bed, from_string=True)
        cds_bedtool = CdsBedToolBuilder().generate_bedtool(gene_db)
        utr3_bedtool = utr3_bedtool.subtract(cds_bedtool).remove_invalid().sort()
        return utr3_bedtool


class StartCodonBedToolBuilder(AbstractBedToolBuilder):
    _feature_name = 'start_codon'

    def _make_bed(self, gene_db):
        f_bed = ""
        for gene_id in gene_db.get_gene_list():
            f_regions = []
            for feature in gene_db.feature_db.children(gene_id, featuretype=self._feature_name):
                feature.stop = feature.start
                f_regions.append(feature)

            merged_f = self._merge_regions(gene_db.feature_db, f_regions)
            renamed_f = self._rename_regions(merged_f, gene_id)
            f_bed += self._create_bed(renamed_f)

        return f_bed


class StopCodonBedToolBuilder(AbstractBedToolBuilder):
    _feature_name = 'stop_codon'

    def _make_bed(self, gene_db):
        f_bed = ""
        for gene_id in gene_db.get_gene_list():
            f_regions = []
            for feature in gene_db.feature_db.children(gene_id, featuretype=self._feature_name):
                feature.start = feature.stop
                feature.stop = feature.stop + 1
                f_regions.append(feature)

            merged_f = self._merge_regions(gene_db.feature_db, f_regions)
            renamed_f = self._rename_regions(merged_f, gene_id)
            f_bed += self._create_bed(renamed_f)

        return f_bed


class BedToolBuilderFactory(object):
    """
    This class act as a factory, which delivers object corresponding to the feature.
    """

    _builders = {
        "gene": GeneBedToolBuilder,
        "exon": ExonBedToolBuilder,
        "intron": IntronBedToolBuilder,
        "CDS": CdsBedToolBuilder,
        "5utr": FivePrimeUtrBedToolBuilder,
        "3utr": ThreePrimeUtrCdsBedToolBuilder,
        "start_codon": StartCodonBedToolBuilder,
        "stop_codon": StopCodonBedToolBuilder
    }

    @staticmethod
    def get_builder(feature):
        return BedToolBuilderFactory._builders[feature]()


