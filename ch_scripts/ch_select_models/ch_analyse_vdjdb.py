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

pairs = vdj.get_alleles(vdjdb, paired_only=False, drop_duplicates=True,
                        change_part_v=[r'\*.*', ''],
                        change_part_mhc=[r'(\*\d{2}).*', r'\1'])
for ch_tmp in ['b', 'a']:
    pairs = ch_funcs.ch_column_value_replace(pairs, 'mhc.{}'.format(ch_tmp), ['HLA-', ''])
pairs['ref.pdb'] = ''


def find_potential_structures(pdb_db, seq_db, paired=True, distinguish_alpha_beta_for_future_drop=True):
    same_v_mhc = pd.DataFrame(columns=seq_db.columns)
    same_v_mhc['ref.pdb'] = ''
    same_v_mhc['pair.type'] = ''
    for index, row in pdb_db.iterrows():
        species = row['species']
        mhc_a = str(row['mhc_a_allele'])
        mhc_b = str(row['mhc_b_allele'])
        if species == 'MusMusculus':
            mhc_a = mouse_mhc_a[mhc_a]
            mhc_b = mouse_mhc_b[mhc_b]
        tcr_v_a = row['tcr_v_alpha_allele']
        tcr_v_b = row['tcr_v_beta_allele']
        antigen = row['antigen_seq']
        cdr3a = row['cdr3a']
        cdr3b = row['cdr3b']
        pdb = row['pdb_id']
        if paired is True:
            chains = copy.deepcopy(seq_db[(seq_db['species'] == species)
                           & (seq_db['antigen.epitope'] == antigen)
                           & (seq_db['mhc.a'] == mhc_a)
                           & (seq_db['mhc.b'] == mhc_b)
                           & (seq_db['v.alpha.segm'] == tcr_v_a)
                           & (seq_db['v.beta.segm'] == tcr_v_b)
                           & ((seq_db['cdr3.alpha'] != cdr3a) | (seq_db['cdr3.beta'] != cdr3b))
                           & (seq_db['cdr3.alpha'].str.len() == len(cdr3a))
                           & (seq_db['cdr3.beta'].str.len() == len(cdr3b))])
            if len(chains) > 0:
                chains.loc[:, ('ref.pdb')] = pdb
                chains.loc[:, ('pair.type')] = 'paired'
                if distinguish_alpha_beta_for_future_drop is True:
                    chains.loc[:, ('alphas')] = '+'
                    chains.loc[:, ('betas')] = '+'
            same_v_mhc = same_v_mhc.append(chains, ignore_index=True)
        else:
            chains_a = copy.deepcopy(seq_db[(seq_db['species'] == species)
                                            & (seq_db['antigen.epitope'] == antigen)
                                            & (seq_db['mhc.a'] == mhc_a)
                                            & (seq_db['mhc.b'] == mhc_b)
                                            & (((seq_db['v.alpha.segm'] == tcr_v_a) & (seq_db['cdr3.alpha'].str.len() == len(cdr3a)) & (seq_db['cdr3.alpha'] != cdr3a))
                                               & ((seq_db['v.beta.segm'] != tcr_v_b) | (seq_db['cdr3.beta'].str.len() != len(cdr3b))))])

            chains_b = copy.deepcopy(seq_db[(seq_db['species'] == species)
                                            & (seq_db['antigen.epitope'] == antigen)
                                            & (seq_db['mhc.a'] == mhc_a)
                                            & (seq_db['mhc.b'] == mhc_b)
                                            & (((seq_db['v.alpha.segm'] != tcr_v_a) | (seq_db['cdr3.alpha'].str.len() != len(cdr3a)))
                                             & ((seq_db['v.beta.segm'] == tcr_v_b) & (seq_db['cdr3.beta'].str.len() == len(cdr3b)) & (seq_db['cdr3.beta'] != cdr3b)))])
            if len(chains_a) > 0:
                chains_a.loc[:, ('ref.pdb')] = pdb
                chains_a.loc[:, ('pair.type')] = 'unpaired.a'
                if distinguish_alpha_beta_for_future_drop is True:
                    chains_a.loc[:, ('alphas')] = '+'
                    chains_a.loc[:, ('betas')] = '-'
            if len(chains_b) > 0:
                chains_b.loc[:, ('ref.pdb')] = pdb
                chains_b.loc[:, ('pair.type')] = 'unpaired.b'
                if distinguish_alpha_beta_for_future_drop is True:
                    chains_b.loc[:, ('alphas')] = '-'
                    chains_b.loc[:, ('betas')] = '+'
            same_v_mhc = same_v_mhc.append(chains_a, ignore_index=True).append(chains_b, ignore_index=True)
    return same_v_mhc

same_v_mhc = find_potential_structures(unique_pdbs, pairs, paired=True).append(find_potential_structures(unique_pdbs, pairs, paired=False), ignore_index=True)
print(len(same_v_mhc))
same_v_mhc.to_csv(os.path.join(BASE_SCRIPTDIR, 'TCR_MHC_antigen_sequences_unpaired_for_modelling.txt'), sep='\t', index=None)

same_v_mhc = same_v_mhc.drop_duplicates(subset=['complex.id', 'alphas'], keep='first').drop_duplicates(subset=['complex.id', 'betas'], keep='first').reset_index(drop=True).drop(columns=['betas', 'alphas'])
same_v_mhc.to_csv(os.path.join(BASE_SCRIPTDIR, 'TCR_MHC_antigen_sequences_unpaired_no_duplicates.txt'), sep='\t', index= None)
print(len(same_v_mhc))

#check antigens
#same_v_mhc['cdr3.beta.len'] = same_v_mhc['cdr3.beta'].str.len()
#same_v_mhc['cdr3.alpha.len'] = same_v_mhc['cdr3.alpha'].str.len()
#group_mhc = same_v_mhc[same_v_mhc['mhc.class']=='MHCI'].groupby(['mhc.a', 'antigen.epitope', 'cdr3.alpha.len', 'cdr3.beta.len', 'pair.type'])['pair.type'].agg(['count']).reset_index()#.value_counts()#.reset_index()#['paired', 'unpaired.alpha', 'unpaired.beta'].sum().reset_index()
#group_mhc_1 = ch_funcs.ch_partial_melt(group_mhc, 'pair.type', change_columns=['count'], save_columns=['mhc.a', 'antigen.epitope', 'cdr3.alpha.len', 'cdr3.beta.len']).fillna(0)#['mhc.a', 'antigen.epitope', 'cdr3.alpha.len', 'cdr3.beta.len']).fillna(0)
#print(group_mhc_1)
#group_mhc.to_csv(os.path.join(BASE_SCRIPTDIR, 'MHCI_info_about_possible_structures.txt'), sep='\t', index=None)


