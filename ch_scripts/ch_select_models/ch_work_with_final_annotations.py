import pandas as pd

from ch_base import *
import ch_db_functions as db_funcs

file = os.path.join(BASE_HOMEDIR, "ch_scripts/modeller/final.annotations.txt")
pdbfiles = os.path.join(BASE_HOMEDIR, "allpdb/")


def get_alleles(input_path, change_part_v=None, change_part_mhc=None):
    inpcols = ['pdb_id', 'mhc_a_allele', 'mhc_b_allele', 'tcr_v_allele', 'antigen_seq', 'mhc_type', 'species', 'tcr_region', 'tcr_region_seq']
    annot = pd.read_table(input_path, usecols=inpcols)[inpcols]
    annot = annot[annot['tcr_region'] == 'CDR3']
    annot = annot.drop_duplicates(subset=['pdb_id', 'tcr_v_allele'], keep='first').reset_index(drop=True)
    annot['chain_type'] = annot['tcr_v_allele'].str[2]

    beta = annot[annot['chain_type']=='B'].drop(columns=['chain_type', 'tcr_region'])
    alpha = annot[annot['chain_type'] == 'A'].drop(columns=['chain_type', 'tcr_region'])

    beta.columns = ['pdb_id', 'mhc_a_allele', 'mhc_b_allele', 'tcr_v_beta_allele', 'antigen_seq', 'mhc_type', 'species', 'cdr3b']
    alpha.columns = ['pdb_id', 'mhc_a_allele', 'mhc_b_allele', 'tcr_v_alpha_allele', 'antigen_seq', 'mhc_type', 'species', 'cdr3a']

    pdbs = alpha.merge(beta, on=['pdb_id', 'mhc_a_allele', 'mhc_b_allele', 'antigen_seq', 'mhc_type', 'species'], how='outer')


    if change_part_v is not None and isinstance(change_part_v, (list, tuple)):
        for ch_tmp in ['beta', 'alpha']:
            pdbs = db_funcs.ch_column_value_replace(pdbs, 'tcr_v_{}_allele'.format(ch_tmp), change_part_v)
    if change_part_mhc is not None and isinstance(change_part_mhc, (list, tuple)):
        for ch_tmp in ['b', 'a']:
            pdbs = db_funcs.ch_column_value_replace(pdbs, 'mhc_{}_allele'.format(ch_tmp), change_part_mhc)

    return pdbs


def get_chains(input_path, var = 1):
    annot = pd.read_table(input_path, usecols=['pdb_id', 'chain_tcr', 'tcr_v_allele', 'chain_mhc_a', 'chain_mhc_b', 'chain_antigen'])
    annot = annot.drop_duplicates(subset=['pdb_id', 'chain_tcr', 'tcr_v_allele'], keep='first').reset_index(drop=True)
    annot['chain_type'] = annot['tcr_v_allele'].str[2]

    pdbs = {pdb: chain.set_index('chain_type').to_dict('index')
            for pdb,
                chain in annot.groupby('pdb_id')}
    if var == 1:
        return pdbs
    else:
        pdbs_var2 = {key: {'chain_alpha':pdbs[key]['A']['chain_tcr'],
                           'chain_beta':pdbs[key]['B']['chain_tcr'],
                           'chain_mhc_a':pdbs[key]['A']['chain_mhc_a'],
                           'chain_mhc_b':pdbs[key]['A']['chain_mhc_b'],
                           'chain_antigen':pdbs[key]['A']['chain_antigen'],
                           'trav_allele':pdbs[key]['A']['tcr_v_allele'],
                           'trbv_allele':pdbs[key]['B']['tcr_v_allele']} for key in pdbs}
    return pdbs_var2


def dicttolist_pdbs(pdbs, type='all', use_chains='all'):
    for pdb in pdbs:
        chains = {}
        for chain in pdbs[pdb]:
            info = pdbs[pdb][chain]
            if use_chains=='all':
                pass
            else:
                if chain != use_chains:
                    continue
            if not pd.isnull(info['chain_tcr']) and (type == 'all' or type == 'only_tcr'):
                chains[info['chain_tcr']] = chain
            if not pd.isnull(info['chain_mhc_b']) and (type == 'all'):
                chains[info['chain_mhc_b']] = 'chain_mhc_b'
            if not pd.isnull(info['chain_mhc_a']) and (type == 'all'):
                chains[info['chain_mhc_a']] = 'chain_mhc_a'
            if not pd.isnull(info['chain_antigen']) and (type == 'all'):
                chains[info['chain_antigen']] = 'chain_antigen'
        if (use_chains=='all' and type == 'all' and len(chains) != 5) or (use_chains!='all' and type == 'all' and len(chains) != 4):
            print(pdb, chains)
            continue
        pdbs[pdb] = chains
    return pdbs


def get_cdr_positions(inpfile = file, region='CDR3'):
    finannot = pd.read_table(inpfile)
    finannot['chain_type'] = finannot['tcr_v_allele'].str[2]
    inpcols = ['pdb_id', 'tcr_region_start', 'tcr_region_end', 'tcr_region_seq', 'chain_type']
    finannot = finannot[finannot['tcr_region'] == region][inpcols].reset_index(drop=True)
    return {pdb: chain.set_index('chain_type').to_dict('index')
     for pdb,
         chain in finannot.groupby('pdb_id')}
