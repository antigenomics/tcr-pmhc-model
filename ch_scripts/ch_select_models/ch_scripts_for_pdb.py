from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#this code was obtained from bioinfopoint.com (http://www.bioinfopoint.com/index.php/code/41-extracting-the-sequence-from-a-protein-structure-using-biopython)
def get_amino_seq_chains(inpfile, use_chains='all'):
    p = PDBParser(PERMISSIVE=1)
    structure = p.get_structure(inpfile, inpfile)
    model = structure[0]
    chains = {}
    for chain in model:
        if chain.id in use_chains or use_chains=='all':
            chains[chain.id] = ""
            for residue in chain:
                if is_aa(residue.get_resname(), standard=True):
                    chains[chain.id] += three_to_one(residue.get_resname())
                else:
                    chains[chain.id] += "X"
    return chains

