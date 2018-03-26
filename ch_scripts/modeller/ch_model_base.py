from modeller import *
from modeller.automodel import *
from modeller import soap_protein_od
from modeller.scripts import complete_pdb

import subprocess
import sys
import os
from shutil import copyfile

def run_blast(filename, dbpath):
    out_file = filename + ".blast"
    cmd = "psiblast -db " + dbpath + " -query " + filename + " -out " + out_file
    # Put stderr and stdout into pipes
    proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    return_code = proc.wait()
    if return_code != 0:
        print("blast subprocess failed with exit code ", return_code)


def get_fasta_info(filename, type):
    inp = open(filename, 'r')
    info = inp.readline().replace('>', '', 1).strip().split('|')
    #print('get_fasta_info', info)
    if type == 'pdb':
        inp.close()
        return(info[-1], info[1])
    else:
        inp.close()
        return(info[-1], info[-2])

def get_pdb(pdb_code, pdb_dir, dir_path):
    # UPDATE PDB repository
    repository_pdb = pdb_dir#"/Users/AlekseyYeliseev/Desktop/AdaptiveImm/ch_scripts/chain_PDBs"
    #divided_pdb_folder = pdb_code[1:3];
    pdb_inp = repository_pdb + "/pdb" + pdb_code + ".ent"
    pdb_out = dir_path + "/pdbs/{0}.pdb".format(pdb_code)
    copyfile(pdb_inp, pdb_out)


def get_float_evalue(string):
    if string.startswith('e'):
        string = string.replace('e', '1e')
    return(float(string))


def modeller_struct_multiple_alignment(aln, output_file):
    for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                        ((1., 0.5, 1., 1., 1., 0.), False, True),
                                        ((1., 1., 1., 1., 1., 0.), True, False)):
        aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50),
                   gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                   alignment_type='tree',  # If 'progresive', the tree is not
                   # computed and all structues will be
                   # aligned sequentially to the first
                   feature_weights=weights,  # For a multiple sequence alignment only
                   # the first feature needs to be non-zero
                   improve_alignment=True, fit=True, write_fit=write_fit,
                   write_whole_pdb=whole, output='ALIGNMENT QUALITY')

    aln.write(file=output_file, alignment_format='PIR')
    return aln


def modeller_alignment_add_sequence(aln, output_file):
    aln.salign(output='', max_gap_length=10,
               gap_function=True,  # to use structure-dependent gap penalty
               alignment_type='PAIRWISE', align_block=len(aln),
               feature_weights=(1., 0., 0., 0., 0., 0.), overhang=99,
               gap_penalties_1d=(-450, 0),
               gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
               similarity_flag=True, local_alignment=True)

    aln.write(file=output_file, alignment_format='PIR')
    return aln


def modeller_search_templates(inpfile, var='default'):
    if var == 'default':
        inp = open(inpfile, 'r')
        template_list = inp.read().strip().split('\n')
        inp.close()
    return template_list


def remove_files(num):
    for n in range(1, num+1):
        file_name1 = "TCR.D0000{0:{fill}{align}4}".format(n, fill=0, align='>')
        file_name2 = "TCR.V9999{0:{fill}{align}4}".format(n, fill=0, align='>')
    os.remove("TCR.sch")
    os.remove("TCR.ini")
    os.remove("TCR.rsr")

