import os
import sys
import pandas as pd
from Bio import PDB
from ch_base import *
from ch_work_with_final_annotations import *
from pdb_splitter import *

FINAL_ANNOTATIONS = os.path.join(BASE_HOMEDIR, "ch_scripts/modeller/final.annotations.txt")
PATH_TO_PDB = os.path.join(BASE_HOMEDIR, "allpdb/")

chaindict = {'A': 'alpha', 'B': 'beta'}

pdbs = get_chains(FINAL_ANNOTATIONS, var=2)


def complex_type_options(type_chains):
    """
    By this function you can manage what chains will be left untouched (other chains will be removed)
    :param type_chains: just write all chains (alpha/beta/mhc)
    :return: dictionary with options
    """
    complex_type = {'chain_alpha': True, 'chain_beta': True, 'chain_mhc_a': True, 'chain_mhc_b': True, 'chain_antigen': True}
    if type_chains == 'all':
        pass
    else:
        if 'alpha' not in type_chains:
            complex_type['chain_alpha']=False
        if 'beta' not in type_chains:
            complex_type['chain_beta']=False
        if 'mhc' not in type_chains:
            complex_type['chain_mhc_a'] = False
            complex_type['chain_mhc_b'] = False
            if 'antigen' not in type_chains:
                complex_type['chain_antigen'] = False
        if 'no_antigen' in type_chains: #just in case..
            complex_type['chain_antigen'] = False

    return complex_type


def rewrite_chain(inpfile, inpchains, outchains):
    pdb_list = []
    chain_chainger = {}
    chain_dict = {}

    for i in range(len(inpchains)):
        chain_chainger[inpchains[i]] = outchains[i]
        chain_dict[outchains[i]] = []
    with open(inpfile, 'r') as inp:
        for line in inp:
            if line.startswith('END'):
                pdb_list.append(line)
                break
            info = list(line)
            info[21] = chain_chainger[info[21]] #chain
            chain_dict[info[21]].append("".join(info)) #atoms in chain

    with open(inpfile, 'w') as out:
        iti = 0
        winfo = 'JUSTRANDOMSTRING'
        for chain in sorted(outchains):
            for line_iti in range(len(chain_dict[chain])):
                workline = chain_dict[chain][line_iti]
                if workline[0:4] == 'ATOM':
                    if winfo != workline[22:26]:
                        iti += 1
                        winfo = workline[22:26]
                    workline_change = list(workline)
                    workline_change[22:26] = '{0:{fill}{align}4}'.format(iti, fill='', align='>')
                    chain_dict[chain][line_iti] = ''.join(workline_change)
            out.writelines(chain_dict[chain])


def get_splitted_pdbs(pdb_path=PATH_TO_PDB, type_chains='all', pdbs=pdbs):
    """
    This script will change structures from pdbs: it will get specific chains from pdb file and rename them:
    A – TCRalpha
    B – TCRbeta
    C – antigen
    D – MHCa
    E – MHCe
    :param pdb_path:
    :param type_chains:
    :param pdbs:
    :return:
    """
    print('Working with {}'.format(type_chains))
    complex_type = complex_type_options(type_chains)
    splitter = PDBSplitter(out_dir='{}_PDBs'.format(type_chains))
    for pdb in pdbs:
        w_chains = ''
        to_change = ''
        w_pdb = pdbs[pdb]
        for chain in sorted(complex_type):
            if complex_type[chain] is True:
                w_chains += w_pdb[chain]
                if chain is 'chain_alpha':
                    to_change += 'A'
                elif chain is 'chain_beta':
                    to_change += 'B'
                elif chain is 'chain_antigen':
                    to_change += 'C'
                elif chain is 'chain_mhc_a':
                    to_change += 'D'
                elif chain is 'chain_mhc_b':
                    to_change += 'E'
        outpdb = splitter.split_by_chains(pdb_path=pdb_path, pdb_id='pdb{}.ent'.format(pdb), chain_letters=w_chains, overwrite=True, outfile_name='pdb{}.ent'.format(pdb))
        rewrite_chain(outpdb, w_chains, to_change)


#for type_chains in ['all', 'tcr.alpha.mhc', 'tcr.beta.mhc']:
#    get_splitted_pdbs(PATH_TO_PDB, type_chains, pdbs)
