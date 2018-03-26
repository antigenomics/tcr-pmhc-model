import re

import pandas as pd

from ch_base import *
import ch_db_functions as db_funcs


vdjdb = os.path.join(BASE_HOMEDIR, "vdjdb-2018-02-03/vdjdb.txt")


def get_alleles(input_path, paired_only=True, drop_duplicates=True, change_part_v=None, change_part_mhc=None):
    inpcols = ['complex.id', 'gene', 'cdr3', 'v.segm', 'j.segm', 'species', 'mhc.a', 'mhc.b', 'mhc.class',
               'antigen.epitope', 'meta']

    annot = pd.read_table(input_path, usecols=inpcols)[inpcols]
    annot = annot[(annot['meta'].str.contains('"structure.id": ""')) | (annot['meta'].str.contains('"structure.id": "NA"'))].drop(columns=['meta'])


    acol = ['complex.id', 'cdr3.alpha', 'v.alpha.segm', 'j.alpha.segm', 'species',
            'mhc.a', 'mhc.b', 'mhc.class', 'antigen.epitope']
    bcol = ['complex.id', 'cdr3.beta', 'v.beta.segm', 'j.beta.segm', 'species',
            'mhc.a', 'mhc.b', 'mhc.class', 'antigen.epitope']

    if not paired_only:
        unpaired = annot[annot['complex.id'] == 0]
        unpaired_alpha = unpaired[unpaired['gene'] == 'TRA'].drop(columns=['gene'])
        unpaired_beta = unpaired[unpaired['gene'] == 'TRB'].drop(columns=['gene'])
        unpaired_alpha.columns = acol
        unpaired_beta.columns = bcol
        unpaired = unpaired_alpha.append(unpaired_beta, ignore_index=True)
        unpaired_ids = [-i for i in range(1, len(unpaired)+1)]
        unpaired['complex.id'] = unpaired_ids
    
    annot = annot[annot['complex.id'] != 0] #remove unpaired TCRs

    alpha = annot[annot['gene'] == 'TRA'].drop(columns=['gene'])
    beta = annot[annot['gene'] == 'TRB'].drop(columns=['gene'])

    alpha.columns = acol
    beta.columns = bcol

    vdjdb = alpha.merge(beta, on=['complex.id', 'species', 'mhc.a', 'mhc.b', 'mhc.class',
                                  'antigen.epitope'],
                        how='outer')

    if not paired_only:
        vdjdb = vdjdb.append(unpaired, ignore_index=True)


    if change_part_v is not None and isinstance(change_part_v, (list, tuple)):
        for ch_tmp in ['beta', 'alpha']:
            vdjdb = db_funcs.ch_column_value_replace(vdjdb, 'v.{}.segm'.format(ch_tmp), change_part_v)
    if change_part_mhc is not None and isinstance(change_part_mhc, (list, tuple)):
        for ch_tmp in ['b', 'a']:
            vdjdb = db_funcs.ch_column_value_replace(vdjdb, 'mhc.{}'.format(ch_tmp), change_part_mhc)

    if drop_duplicates:
        vdjdb = vdjdb.drop_duplicates(subset=['cdr3.alpha', 'v.alpha.segm', 'j.alpha.segm',
                                              'cdr3.beta', 'v.beta.segm', 'j.beta.segm',
                                              'species', 'mhc.a', 'mhc.b',
                                              'antigen.epitope'],
                                      keep='first')
    return vdjdb
