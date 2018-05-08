import os
import sys
import Bio.PDB
from ch_base import *


class Found(Exception): pass


class NotDisordered(Bio.PDB.Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == 'A'


pdb_parser = Bio.PDB.PDBParser(QUIET = True)


def get_structure_and_model(inppath, type):
    structure = pdb_parser.get_structure(type, inppath)
    model = structure[0]
    return structure, model


def get_atoms_from_chain(model, chain):
    atoms = []
    for res in model[chain]:
        if res.get_id()[0] == ' ':
            atoms.append(res['CA'])
    return atoms


def select_atoms(atoms, inprange):
    return atoms[inprange[0]:inprange[1]]


def select_surrounding_atoms(atoms, inprange, depth=3):
    downrange, uprange = (inprange[0]-depth, inprange[0]), (inprange[1]+1, inprange[1]+1+depth)
    return select_atoms(atoms, downrange)+select_atoms(atoms, uprange)


#made in past, divided to get_structure and get_atoms_from_chain
def ch_get_structure(inppath, type, chain):
    structure = pdb_parser.get_structure(type, inppath)
    model = structure[0]
    atoms = get_atoms_from_chain(model, chain)
    return atoms, model, structure


def superimpose_two(ref_atoms, sample_atoms, sample_model, sample_structure, outfile):
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure)
    io.save(outfile, select=NotDisordered())


def superimpose_sh(ref_pdb, sample_pdb, ref_chain, sample_chain, inppath, outpath, antigen):
    ref_atoms, _, _ = ch_get_structure(os.path.join(inppath, ref_pdb), 'reference', ref_chain)
    sample_atoms, sample_model, sample_structure = ch_get_structure(
        os.path.join(inppath, sample_pdb), 'sample', sample_chain)

    outfile = os.path.join(outpath, '{}_align_to_{}.pdb'.format(sample_pdb, ref_pdb))
    superimpose_two(ref_atoms, sample_atoms, sample_model, sample_structure, outfile)


def superimpose_(pdbdict, outpath, overwrite=True):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    for antigen in pdbdict:
        for idpdb1 in range(len(pdbdict[antigen])):
            for idpdb2 in range(idpdb1+1, len((pdbdict[antigen]))):
                pdb1 = sorted(pdbdict[antigen])[idpdb1]
                pdb2 = sorted(pdbdict[antigen])[idpdb2]
                if overwrite == False:
                    if os.path.exists(os.path.join(outpath, antigen, pdb2+"_align_"+pdb1+".pdb")):
                        continue
                for chain in pdbdict[antigen][pdb1]:
                    chain1 = pdbdict[antigen][pdb1][chain]['pdb_antigen']
                    break
                for chain in pdbdict[antigen][pdb2]:
                    chain2 = pdbdict[antigen][pdb2][chain]['pdb_antigen']
                    break


#refstructure, refmodel = get_structure_and_model(os.path.join(BASE_HOMEDIR, 'ch_scripts/ch_select_models/all_PDBs/pdb1ao7.ent'), 'reference')
#samplestructure, samplemodel = get_structure_and_model(os.path.join(BASE_HOMEDIR, 'ch_scripts/ch_select_models/all_PDBs/pdb3h9s.ent'), 'sample')
#refatoms_a, refatoms_b = get_atoms_from_chain(refmodel, 'A'), get_atoms_from_chain(refmodel, 'B')
#refatoms = select_surrounding_atoms(refatoms_a, (87, 100), 10) + select_surrounding_atoms(refatoms_b, (88, 104), 10)

#sampleatoms_a, sampleatoms_b = get_atoms_from_chain(samplemodel, 'A'), get_atoms_from_chain(samplemodel, 'B')
#sampleatoms = select_surrounding_atoms(sampleatoms_a, (87, 100), 10) + select_surrounding_atoms(sampleatoms_b, (90, 106), 10)

#superimpose_two(refatoms, sampleatoms, samplemodel, samplestructure, os.path.join(BASE_HOMEDIR, 'ch_scripts/ch_alignstruct/pdb3h9s_to_1ao7.ent'))