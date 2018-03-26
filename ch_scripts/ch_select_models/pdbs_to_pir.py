import os
from ch_base import *
from ch_work_with_final_annotations import *
import ch_scripts_for_pdb as ch_pdb_funcs

file = os.path.join(BASE_HOMEDIR, "ch_scripts/modeller/final.annotations.txt")

pdbfiles = os.path.join(BASE_SCRIPTDIR, "all_PDBs/")

pdbs = get_chains(file, var=2)

with open('all_pdb_sequences.txt', 'w') as out: #A/B - alpha/beta; C – antigen; D – mhc_a; E – mhc.b
    out.write('pdb\tchain\tseq\n')
    for pdb in pdbs:
        pdb_chains = ch_pdb_funcs.get_amino_seq_chains(os.path.join(pdbfiles, 'pdb{}.ent'.format(pdb)))
        for chain in sorted(pdb_chains):
            out.write('{}\t{}\t{}\n'.format(pdb, chain, pdb_chains[chain].replace('X', '')))

