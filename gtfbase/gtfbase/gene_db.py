"""
This module will contain methods to create a Gene Dictionary (of type
DefaultOrderedDictionary) and a FeatureDB object from a GTF file.
"""
import logging

import gffutils
from default_ordered_dictionary import DefaultOrderedDict


class GeneDB(object):

    def __init__(self, gtf_path):
        """
        Constructor for GeneDB.
        """
        self._gtf_path = gtf_path
        self._feature_db = None
        self._gene_dict = None

    @property
    def feature_db(self):
        if self._feature_db is None:
            self._feature_db = gffutils.create_db(self._gtf_path, dbfn=":memory:",
                                                  merge_strategy='merge',
                                                  force=True,
                                                  disable_infer_transcripts=True,
                                                  disable_infer_genes=True,
                                                  keep_order=True)
        return self._feature_db

    @property
    def gene_dict(self):
        if self._gene_dict is None:
            self._gene_dict = self._create_gene_dict()
        return self._gene_dict

    def features_available(self):
        return list(self.feature_db.featuretypes())

    def _create_gene_dict(self):
        '''
        Store each feature line db.all_features() as a dict of dicts
        '''
        gene_dict = DefaultOrderedDict(lambda: DefaultOrderedDict(lambda: DefaultOrderedDict(list)))
        for line_no, feature in enumerate(self.feature_db.all_features()):
            gene_ids = feature.attributes['gene_id']
            feature_type = feature.featuretype
            if feature_type == 'gene':
                if len(gene_ids)!=1:
                    logging.warning('Found multiple gene_ids on line {} in gtf'.format(line_no))
                    break
                else:
                    gene_id = gene_ids[0]
                    gene_dict[gene_id]['gene'] = feature
            else:
                transcript_ids = feature.attributes['transcript_id']

                for gene_id in gene_ids:
                    for transcript_id in transcript_ids:
                        gene_dict[gene_id][transcript_id][feature_type].append(feature)
        return gene_dict

    def get_gene_list(self):
        return list(set(self.gene_dict.keys()))