import re
import copy
from functools import reduce

import pandas as pd

from ch_base import *
import ch_work_with_final_annotations as annot
import ch_work_with_vdjdb as vdj
import ch_db_functions as ch_funcs

final_annotations = os.path.join(BASE_HOMEDIR, "ch_scripts/modeller/final.annotations.txt")
pdbfiles = os.path.join(BASE_HOMEDIR, "allpdb/")
vdjdb = os.path.join(BASE_HOMEDIR, "vdjdb-2018-02-03/vdjdb.txt")


pdbs = annot.get_alleles(final_annotations,
                         change_part_v=[r'\*.*', ''],
                         change_part_mhc=[r'(\*\d{2}).*', r'\1'])

pdbs['species'] = pdbs['species'].map({'Homo_sapiens':'HomoSapiens', 'Mus_musculus':'MusMusculus'})


unique_pdbs = pdbs.drop_duplicates(subset=['mhc_a_allele', 'mhc_b_allele', 'tcr_v_beta_allele', 'tcr_v_alpha_allele',
                                           'antigen_seq', 'mhc_type', 'species', 'cdr3b', 'cdr3a'],
                                   keep='first').reset_index(drop=True)

pairs = vdj.get_alleles(vdjdb, paired_only=True, drop_duplicates=True,
                        change_part_v=[r'\*.*', ''],
                        change_part_mhc=[r'(\*\d{2}).*', r'\1'])
for ch_tmp in ['b', 'a']:
    pairs = ch_funcs.ch_column_value_replace(pairs, 'mhc.{}'.format(ch_tmp), ['HLA-', ''])
pairs['ref.pdb'] = ''

pdbs_epitopes = set(unique_pdbs['antigen_seq'])
print(len(unique_pdbs), sorted(list(pdbs_epitopes)))
pairs_epitopes = set(pairs[pairs['mhc.class']=='MHCI']['antigen.epitope'])
print(len(pdbs_epitopes & pairs_epitopes), len(pairs_epitopes), len(pdbs_epitopes))