import pandas as pd
import csv
from ch_base import *
from ch_work_with_final_annotations import get_cdr_positions, get_chains
import ch_db_functions as db_funcs
import ch_scripts_for_pdb as ch_pdb_funcs

file = os.path.join(BASE_HOMEDIR, "ch_scripts/modeller/final.annotations.txt")
model_files = os.path.join(BASE_SCRIPTDIR, "TCR_MHC_antigen_sequences_unpaired_no_duplicates.txt")
pdbfiles = 'all_pdb_sequences.txt'


#What we need for paired: Combine V J and get cdr3! Yay :3
#What we need for unpaired: THE SAME!


cdr1s = get_cdr_positions(file, 'CDR1')
cdr2s = get_cdr_positions(file, 'CDR2')
cdr3s = get_cdr_positions(file, 'CDR3')
cdrs = [cdr1s, cdr2s, cdr3s]


pdbs = {}
with open(pdbfiles, 'r') as inp:
    reader = csv.DictReader(inp, delimiter='\t')
    for row in reader:
        if row['pdb'] not in pdbs:
            pdbs[row['pdb']] = {}
        pdbs[row['pdb']][row['chain']] = row['seq']

def create_seq(pdb, cdr3, chain):
    pdb_sequence = pdbs[pdb][chain].replace('X', '')
    start = cdr3s[pdb][chain]['tcr_region_start']
    end = cdr3s[pdb][chain]['tcr_region_end']
    pdb_sequence='{}{}{}'.format(pdb_sequence[:start], cdr3, pdb_sequence[start+len(cdr3):])
    return pdb_sequence



def cdr_pos(inp_cdr):
    if inp_cdr is '':
        return ''
    else:
        return '{}_{}'.format(inp_cdr['tcr_region_start'], inp_cdr['tcr_region_end'])

def write_sequences_pir(inpfile=model_files, w_type='paired', w_species='HomoSapiens', outfile='paired_sequences_for_modeller_homosapiens.txt'):
    with open(outfile, 'w') as out:
        with open(inpfile, 'r') as inp:
            reader = csv.DictReader(inp, delimiter='\t')
            for row in reader:
                if row['species'] != w_species:
                    continue
                ref_pdb = row['ref.pdb']
                seq_c = pdbs[ref_pdb]['C']
                seq_d = pdbs[ref_pdb]['D']
                seq_e = pdbs[ref_pdb]['E']

                if row['pair.type'] == 'paired':
                    cdr1a, cdr2a, cdr3a = (a[ref_pdb]['A'] for a in cdrs)
                    cdr1b, cdr2b, cdr3b = (a[ref_pdb]['B'] for a in cdrs)
                elif row['pair.type'] == w_type:
                    if w_type == 'unpaired.a':
                        cdr1a, cdr2a, cdr3a = (a[ref_pdb]['A'] for a in cdrs)
                        cdr1b, cdr2b, cdr3b = '', '', ''

                    elif w_type == 'unpaired.b':
                        cdr1a, cdr2a, cdr3a = '', '', ''
                        cdr1b, cdr2b, cdr3b = (a[ref_pdb]['B'] for a in cdrs)
                else:
                    continue

                if w_type == 'paired' or w_type == 'unpaired.a':
                    seq_a = create_seq(ref_pdb, row['cdr3.alpha'], 'A') + '/'
                else:
                    seq_a = ''
                if w_type == 'paired' or w_type == 'unpaired.b':
                    seq_b = create_seq(ref_pdb, row['cdr3.beta'], 'B') + '/'
                else:
                    seq_b = ''

                out.write('>{}|{}|{}|{}|{}|{}|{}|{}|\n{}\n'.format(row['complex.id'], ref_pdb,
                                                                   cdr_pos(cdr1a), cdr_pos(cdr1b),
                                                                   cdr_pos(cdr2a), cdr_pos(cdr2b),
                                                                   cdr_pos(cdr3a), cdr_pos(cdr3b),
                                                       'sequence:TCR:::::::0.00: 0.00'))
                out.write('{}{}{}/{}/{}\n'.format(seq_a, seq_b, seq_c, seq_d, seq_e))


for species in ['HomoSapiens', 'MusMusculus']:
    outfile = 'paired_sequences_for_modeller_{}.txt'.format(species.lower())
    outfile_a = 'alpha_sequences_for_modeller_{}.txt'.format(species.lower())
    outfile_b = 'beta_sequences_for_modeller_{}.txt'.format(species.lower())
    for w_type, w_outfile in (['paired', outfile], ['unpaired.a', outfile_a], ['unpaired.b', outfile_b]):
      write_sequences_pir(w_type=w_type, w_species=species, outfile=w_outfile)