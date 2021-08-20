"""
This module will contain methods to generate a BedTool object from a given
Gene Dictionary.
"""
import pybedtools
import re
import numpy as num
import functools


class AbstractBedToolBuilder(object):
    """
    Act as a Parent class for other bedtool classes.
    """
    _feature_name = None

    def _make_bed(self, gene_db):
        """
        Helper function for generate_bedtool method.
        :param gene_db: geneDB instance
        :return: The bed file content string, will be used to generate a bedtool.
        :rtype: str
        """

        # Todo: if self._feature_name is None: error dena ha

        f_bed = ""
        for transcript_id in gene_db.get_gene_list():
            f_regions = []
            for feature in gene_db.gene_dict[transcript_id].keys():
                if feature == 'gene':
                    continue
                f = list(gene_db.gene_dict[transcript_id][feature][self._feature_name])
                f_regions += f

            merged_f = self._merge_regions(gene_db.feature_db, f_regions)
            renamed_f = self._rename_regions(merged_f, transcript_id)
            f_bed += self._create_bed(renamed_f)

        return f_bed

    @functools.lru_cache(maxsize=5)
    def generate_bedtool(self, gene_db):
        """
        Provides the desired feature bedtool instance.

        :param gene_db: geneDB instance.
        :return: bedtool instance.
        """
        return self._generate_bedtool_helper(gene_db)

    def _generate_bedtool_helper(self, gene_db):
        """
        Helper function for generate_bedtool()

        :param gene_db: geneDB instance.
        :return: bedtool instance.
        """
        f_bed = self._make_bed(gene_db)
        f_bedtool = pybedtools.BedTool(f_bed, from_string=True)
        f_bedtool = f_bedtool.remove_invalid().sort()
        f_bedtool = self._add_score_field(f_bedtool)
        return f_bedtool

    def _add_score_field(self, f_bedtool):
        """

        :param f_bedtool:
        :return:
        """
        cnt = 0
        for data in f_bedtool.features():
            cnt += 1
            data[4] = cnt
        return f_bedtool

    def _create_bed(self, regions, bedtype='0'):
        '''
        Create bed from list of regions
        :param bedtype: 0 or 1
            0-Based or 1-based coordinate of the BED
        :param regions: Regions for which the bed will be created
        :return: bed string
        '''
        bedstr = ''
        for region in regions:
            assert len(region.attributes['transcript_id']) == 1
            ## GTF start is 1-based, so shift by one while writing
            ## to 0-based BED format
            if bedtype == '0':
                start = region.start - 1
            else:
                start = region.start

            bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(region.chrom,
                                                        start,
                                                        region.stop,
                                                        # re.sub('\.\d+', '', region.attributes['transcript_id'][0]),
                                                        # region.attributes['transcript_id'][0]+'_'+region.attributes['gene_id'][0]+'_'+region.attributes['gene_name'][0],
                                                        region.attributes['transcript_id'][0],
                                                        # + '_' +
                                                        # region.attributes['gene_id'][0],
                                                        # + '_' +
                                                        # region.attributes['gene_name'][0],
                                                        '.',
                                                        region.strand)
        return bedstr

    def _rename_regions(self, regions, transcript_id):
        regions = list(regions)
        if len(regions) == 0:
            return []
        for region in regions:
            region.attributes['transcript_id'] = transcript_id
        return regions

    def _merge_regions(self, db, regions):
        """
        This methood is used to merge regions.
        """
        if len(regions) == 0:
            return []
        merged = db.merge(sorted(list(regions), key=lambda x: x.start))
        return merged

    def _merge_regions_nostrand(self, db, regions):
        """
        This methood is used to merge regions while ignoring the strand.
        """
        if len(regions) == 0:
            return []
        merged = db.merge(sorted(list(regions), key=lambda x: x.start), ignore_strand=True)
        return merged


class GeneBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build the Bedtool instance for gene feature.
    """
    _feature_name = 'gene'

    def _make_bed(self, gene_db):
        gene_list = []
        for transcript_id in gene_db.get_gene_list():
            gene_list.append(gene_db.gene_dict[transcript_id]['gene'])

        gene_bed = self._create_bed(gene_list)

        return gene_bed

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


class CdsBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build the Bedtool instance for cds feature.
    """
    _feature_name = 'CDS'


class ExonBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build the Bedtool instance for exon feature.
    """
    _feature_name = 'exon'


class IntronBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build the Bedtool instance for intron feature.
    """
    _feature_name = 'exon'

    def _make_bed(self, gene_db):
        intron_bed = ''

        for transcript_id in gene_db.get_gene_list():
            f_regions = []
            intron_regions = []
            for feature in gene_db.gene_dict[transcript_id].keys():
                if feature == 'gene':
                    continue
                f = list(gene_db.gene_dict[transcript_id][feature][self._feature_name])
                merged_f = self._merge_regions(gene_db.feature_db, f)
                introns = gene_db.feature_db.interfeatures(merged_f)
                f_regions += f
                intron_regions += introns

            merged_introns = self._merge_regions(gene_db.feature_db, intron_regions)
            renamed_introns = self._rename_regions(merged_introns, transcript_id)

            intron_bed += self._create_bed(renamed_introns)

        return intron_bed


class FivePrimeUtrBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build the Bedtool instance for 5'utr feature.
    """
    _feature_name = 'five_prime_utr'

    def _generate_bedtool_helper(self, gene_db):
        utr5_bed = self._make_bed(gene_db)
        utr5_bedtool = pybedtools.BedTool(utr5_bed, from_string=True)
        cds_bedtool = CdsBedToolBuilder()._generate_bedtool_helper(gene_db)
        utr5_bedtool = utr5_bedtool.subtract(cds_bedtool).remove_invalid().sort()
        return utr5_bedtool


class ThreePrimeUtrCdsBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build the Bedtool instance for 3'utr feature.
    """
    _feature_name = 'three_prime_utr'

    def _generate_bedtool_helper(self, gene_db):
        utr3_bed = self._make_bed(gene_db)
        utr3_bedtool = pybedtools.BedTool(utr3_bed, from_string=True)
        cds_bedtool = CdsBedToolBuilder()._generate_bedtool_helper(gene_db)
        utr3_bedtool = utr3_bedtool.subtract(cds_bedtool).remove_invalid().sort()
        return utr3_bedtool


class StartCodonBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build the Bedtool instance for start codon feature.
    """
    _feature_name = 'start_codon'

    def _make_bed(self, gene_db):
        f_bed = ""
        for transcript_id in gene_db.get_gene_list():
            f_regions = []
            for feature in gene_db.feature_db.children(transcript_id, featuretype=self._feature_name):
                feature.stop = feature.start
                f_regions.append(feature)

            merged_f = self._merge_regions(gene_db.feature_db, f_regions)
            renamed_f = self._rename_regions(merged_f, transcript_id)
            f_bed += self._create_bed(renamed_f)

        return f_bed


class StopCodonBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build the Bedtool instance for stop codon feature.
    """
    _feature_name = 'stop_codon'

    def _make_bed(self, gene_db):
        f_bed = ""
        for transcript_id in gene_db.get_gene_list():
            f_regions = []
            for feature in gene_db.feature_db.children(transcript_id, featuretype=self._feature_name):
                feature.start = feature.stop
                feature.stop = feature.stop + 1
                f_regions.append(feature)

            merged_f = self._merge_regions(gene_db.feature_db, f_regions)
            renamed_f = self._rename_regions(merged_f, transcript_id)
            f_bed += self._create_bed(renamed_f)

        return f_bed


class ThreeUtrExonBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build a Bedtool instance for all the exons present in 3'utr.
    """

    def _generate_bedtool_helper(self, gene_db):
        '''
        Uses 3'utr and exon bedtool to generate bedtool for all the exons present in 3'utr.
        '''
        utr3_bedtool = BedToolBuilderFactory.get_builder("3utr")._generate_bedtool_helper(gene_db)
        exon_bedtool = BedToolBuilderFactory.get_builder("exon")._generate_bedtool_helper(gene_db)
        exon_df = exon_bedtool.to_dataframe()
        cnt = 0
        bedstr = ""
        exon_list = exon_df["start"]
        # print(exon_list)
        for utr3 in utr3_bedtool.features():
            start = int(utr3[1])
            stop = int(utr3[2])
            # first  = bisect.bisect_right(exon_list, start, lo=0, hi=len(exon_list))
            first = num.searchsorted(exon_list, start, side='left')
            #   e_start =start
            # Todo: ask what we have to do if an exon has already started but utr is starting.
            for i in range(first, len(exon_list)):
                e_start = exon_df["start"][i]
                e_stop = exon_df["end"][i]

                if e_stop > stop:
                    break
                cnt = cnt + 1
                bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(exon_df["chrom"][i],
                                                            e_start,
                                                            e_stop,
                                                            exon_df["name"][i],
                                                            cnt,
                                                            exon_df["strand"][i])

            if e_start < stop:
                cnt = cnt + 1
                e_stop = stop + 3  # exons are cut at the end of
                bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(exon_df["chrom"][i],
                                                            e_start,
                                                            e_stop,
                                                            exon_df["name"][i],
                                                            cnt,
                                                            exon_df["strand"][i])

        f_bedtool = pybedtools.BedTool(bedstr, from_string=True)
        return f_bedtool


class ThreeUtrIntronBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build a Bedtool instance for all the introns present in 3'utr.
    """

    def _generate_bedtool_helper(self, gene_db):
        '''
        Uses 3'utr and intron bedtool to generate bedtool for all the introns present in 3'utr.
        '''
        utr3_bedtool = BedToolBuilderFactory.get_builder("3utr")._generate_bedtool_helper(gene_db)
        intron_bedtool = BedToolBuilderFactory.get_builder("intron")._generate_bedtool_helper(gene_db)
        intron_df = intron_bedtool.to_dataframe()
        cnt = 0
        bedstr = ""
        intron_list = intron_df["start"]
        # print(exon_list)
        for utr3 in utr3_bedtool.features():
            start = int(utr3[1])
            stop = int(utr3[2])
            # first  = bisect.bisect_right(exon_list, start, lo=0, hi=len(exon_list))
            first = num.searchsorted(intron_list, start, side='left')
            #   e_start =start
            # Todo: ask what we have to do if an exon has already started but utr is starting.
            for i in range(first, len(intron_list)):
                i_start = intron_df["start"][i]
                i_stop = intron_df["end"][i]

                if i_stop > stop:
                    break
                cnt = cnt + 1
                bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(intron_df["chrom"][i],
                                                            i_start,
                                                            i_stop,
                                                            intron_df["name"][i],
                                                            cnt,
                                                            intron_df["strand"][i])

        if i_start < stop:
            i_stop = stop + 3  # exons are cut at the end of
            cnt = cnt + 1
            bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(intron_df["chrom"][i],
                                                        i_start,
                                                        i_stop,
                                                        intron_df["name"][i],
                                                        cnt,
                                                        intron_df["strand"][i])

        f_bedtool = pybedtools.BedTool(bedstr, from_string=True)
        return f_bedtool


class FiveUtrExonBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build a Bedtool instance for all the exons present in 5'utr.
    """

    def _generate_bedtool_helper(self, gene_db):
        '''
        Uses 5'utr and exon bedtool to generate bedtool for all the exons present in 5'utr.
        '''
        utr5_bedtool = BedToolBuilderFactory.get_builder("5utr")._generate_bedtool_helper(gene_db)
        exon_bedtool = BedToolBuilderFactory.get_builder("exon")._generate_bedtool_helper(gene_db)
        exon_df = exon_bedtool.to_dataframe()
        cnt = 0
        bedstr = ""
        exon_list = exon_df["start"]
        # print(exon_list)
        for utr5 in utr5_bedtool.features():
            start = int(utr5[1])
            stop = int(utr5[2])
            # first  = bisect.bisect_right(exon_list, start, lo=0, hi=len(exon_list))
            first = num.searchsorted(exon_list, start, side='left')
            #   e_start =start
            # Todo: ask what we have to do if an exon has already started but utr is starting.
            for i in range(first, len(exon_list)):
                e_start = exon_df["start"][i]
                e_stop = exon_df["end"][i]

                if e_stop > stop:
                    break
                cnt = cnt + 1
                bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(exon_df["chrom"][i],
                                                            e_start,
                                                            e_stop,
                                                            exon_df["name"][i],
                                                            cnt,
                                                            exon_df["strand"][i])

            if e_start < stop:
                e_stop = stop + 3  # exons are cut at the end of
                cnt = cnt + 1
                bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(exon_df["chrom"][i],
                                                            e_start,
                                                            e_stop,
                                                            exon_df["name"][i],
                                                            cnt,
                                                            exon_df["strand"][i])

        f_bedtool = pybedtools.BedTool(bedstr, from_string=True)
        return f_bedtool


class FiveUtrIntronBedToolBuilder(AbstractBedToolBuilder):
    """
    This Class will build a Bedtool instance for all the introns present in 5'utr.
    """

    def _generate_bedtool_helper(self, gene_db):
        '''
        Uses 5'utr and intron bedtool to generate bedtool for all the introns present in 5'utr.
        '''
        utr5_bedtool = BedToolBuilderFactory.get_builder("5utr")._generate_bedtool_helper(gene_db)
        intron_bedtool = BedToolBuilderFactory.get_builder("intron")._generate_bedtool_helper(gene_db)
        intron_df = intron_bedtool.to_dataframe()
        cnt = 0
        bedstr = ""
        intron_list = intron_df["start"]
        # print(exon_list)
        for utr5 in utr5_bedtool.features():
            start = int(utr5[1])
            stop = int(utr5[2])
            first = num.searchsorted(intron_list, start, side='left')
            #   e_start =start
            # Todo: ask what we have to do if an exon has already started but utr is starting.
            for i in range(first, len(intron_list)):
                i_start = intron_df["start"][i]
                i_stop = intron_df["end"][i]

                if i_stop > stop:
                    break
                cnt = cnt + 1
                bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(intron_df["chrom"][i],
                                                            i_start,
                                                            i_stop,
                                                            intron_df["name"][i],
                                                            cnt,
                                                            intron_df["strand"][i])

            if i_start < stop:
                i_stop = stop + 3  # exons are cut at the end of
                cnt = cnt + 1
                bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(intron_df["chrom"][i],
                                                            i_start,
                                                            i_stop,
                                                            intron_df["name"][i],
                                                            cnt,
                                                            intron_df["strand"][i])

        f_bedtool = pybedtools.BedTool(bedstr, from_string=True)
        return f_bedtool


class BedToolBuilderFactory(object):
    """
    The class act as a factory, which delivers object corresponding to the feature.
    """

    _builders = {
        "gene": GeneBedToolBuilder,
        "exon": ExonBedToolBuilder,
        "intron": IntronBedToolBuilder,
        "CDS": CdsBedToolBuilder,
        "5utr": FivePrimeUtrBedToolBuilder,
        "3utr": ThreePrimeUtrCdsBedToolBuilder,
        "start_codon": StartCodonBedToolBuilder,
        "stop_codon": StopCodonBedToolBuilder,
        "3utr_exon": ThreeUtrExonBedToolBuilder,
        "3utr_intron": ThreeUtrIntronBedToolBuilder,
        "5utr_exon": FiveUtrExonBedToolBuilder,
        "5utr_intron": FiveUtrIntronBedToolBuilder
    }

    @staticmethod
    def get_builder(feature):
        """
        The method provides the desired bedtoolfactory object.

        :param feature: The feature for which the factory is required.
        :type feature: str
        :return: BedtoolFactory instance.
        """
        return BedToolBuilderFactory._builders[feature]()
