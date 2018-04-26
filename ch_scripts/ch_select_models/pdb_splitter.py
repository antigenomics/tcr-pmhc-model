from Bio import PDB
from Bio.PDB.PDBIO import Select
from ch_base import *

# code from stackoverflow

class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


class RegionSelect(Select):
    def __init__(self, chain, region_start, region_end):
        self.chain = chain
        self.region_start = region_start
        self.region_end = region_end

    def accept_chain(self, inp_chain):
        if inp_chain.get_id() == self.chain:
            return 1
        else:
            return 0

    def accept_residue(self, residue): #I can rewrite PDBIO.save() function so that it will take residues by
        #  their position in list, not by their resseq number.)
        if residue.get_id()[1] in range(self.region_start, self.region_end):
            return 1
        else:
            return 0


class PDBSplitter:
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

    def split_by_chains(self, pdb_path, pdb_id, chain_letters, overwrite=False, struct=None, outfile_name='pdb.outfile.ent'):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: path to the crystal structure
        :param pdb_id: pdb_id
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """

        chain_letters = [chain.upper() for chain in chain_letters]

        # Input/output files
        out_name = outfile_name
        out_path = os.path.join(self.out_dir, out_name)


        print("OUT PATH:", out_path)
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain{} {} of '{}' already extracted to '{}'.".format(plural,
                                                                         ", ".join(chain_letters),
                                                                         pdb_id, out_name))
            return out_path

        print("Extracting chain{} {} from {}...".format(plural,
                                                        ", ".join(chain_letters),
                                                        pdb_id))

        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, os.path.join(pdb_path, pdb_id))
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters))

        return out_path

    def split_by_region(self, pdb_path, pdb_id, chain, region_start, region_end, overwrite=False, struct=None, outfile_name='pdb.outfile.ent'):
        """

        :param pdb_path: path to the crystal structure
        :param pdb_id: pdb_id
        :param chain: chain id
        :param region_start:
        :param region_end:
        :param overwrite:
        :param struct:
        :param outfile_name:
        :return:
        """
        # Input/output files
        out_name = outfile_name
        out_path = os.path.join(self.out_dir, out_name)

        print("OUT PATH:", out_path)
        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Structure '{}' has been already extracted to '{}'.".format(pdb_id, out_path))
            return out_path
        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, os.path.join(pdb_path, pdb_id))
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=RegionSelect(chain_letter, region_start, region_end))


def split_tcr(pdb_id, tcr_info, pdbdir, outdir):
    """
     This script will split tcr by regions written in tcr_info
    :param pdb_id: name of pdb file
    :param tcr_info: dictionary with structure {chain:{region:[region_start, region_end]}}
    :param pdbdir: path to input pdb
    :param outdir: path to output files
    :return:
    """
    splitter = PDBSplitter(out_dir=outdir)
    for chain in tcr_info:
        for region in tcr_info[chain]:
            region_start = int(tcr_info[chain][region][0])
            region_end = int(tcr_info[chain][region][1])
            splitter.split_by_region(pdb_path=pdbdir, pdb_id=pdb_id, chain=chain,
                                     region_start=region_start, region_end=region_end,
                                     overwrite=True, outfile_name='{}_{}_{}.ent'.format(pdb_id, chain, region))

