import os
import sys
import pandas as pd
from Bio import PDB
from ch_base import *
from ch_work_with_final_annotations import *

file = os.path.join(BASE_HOMEDIR, "ch_scripts/modeller/final.annotations.txt")
pdbfiles = os.path.join(BASE_HOMEDIR, "allpdb/")


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


# code from stackoverflow

class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
            crdir("chain_PDBs")
        else:
            crdir(out_dir)
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, pdb_id, chain_letters, overwrite=False, struct=None):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: path to the crystal structure
        :param pdb_id: pdb_id
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """
        #print(chain_letters)
        chain_letters = [chain.upper() for chain in chain_letters]

        # Input/output files
        full_pdb_id = 'pdb{}.ent'.format(pdb_id)
        out_name = "pdb{}.ent".format(pdb_id)
        out_path = os.path.join(self.out_dir, out_name)


        print("OUT PATH:", out_path)
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain%s %s of '%s' already extracted to '%s'." %
                  (plural, ", ".join(chain_letters), full_pdb_id, out_name))
            return out_path

        print("Extracting chain%s %s from %s..." % (plural,
                                                    ", ".join(chain_letters), pdb_id))

        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(full_pdb_id, os.path.join(pdb_path, full_pdb_id))
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters))

        return out_path


chaindict = {'A':'alpha', 'B':'beta'}

pdbs = get_chains(file, var=2)


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



def get_splitted_pdbs(pdb_path=pdbfiles, type_chains='all', pdbs=pdbs):
    print('Working with {}'.format(type_chains))
    complex_type = complex_type_options(type_chains)
    splitter = ChainSplitter(out_dir='{}_PDBs'.format(type_chains))
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
        outpdb = splitter.make_pdb(pdb_path=pdb_path, pdb_id=pdb, chain_letters=w_chains, overwrite=True)
        rewrite_chain(outpdb, w_chains, to_change)


for type_chains in ['all', 'tcr.alpha.mhc', 'tcr.beta.mhc']:
    get_splitted_pdbs(pdbfiles, type_chains, pdbs)