import os

BASE_SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))
BASE_HOMEDIR = BASE_SCRIPTDIR[:BASE_SCRIPTDIR.rfind('AdaptiveImm')+len('AdaptiveImm')+1]


organism = 'homosapiens'
chain = 'paired'
if chain == 'paired' or chain == 'all':
    inp_pdb_dir = 'all_PDBs'
else:
    inp_pdb_dir = 'tcr.{}.mhc_PDBs'.format(chain)
inp_pdb_dir = os.path.join(BASE_HOMEDIR, 'ch_scripts/ch_select_models', inp_pdb_dir)

mode = "test"#"test_pdb_and_indiscr_antigen"#"indiscr_antigen"#"test_pdb_and_indiscr_antigen"#"indiscr_antigen"#'normal';'test_pdb'

if "test" in mode:
    to_add = 'test'
else:
    to_add = 'real'

workdir = "{}_{}_sequences_{}_seqs_new_var".format(organism, chain, to_add)

num_proc = 2

num_models = 10

sequences = os.path.join(BASE_HOMEDIR, 'ch_scripts/ch_select_models/{}_sequences_for_modeller_{}.txt'.format(chain, organism))#"~/AdaptiveImm/imgt_work/fasta_for_modeller/vdjdb_musmusculus_beta.fasta"#"~/Desktop/AdaptiveImm/imgt_work/fasta_for_modeller/vdjdb_{}_pdb_{}.fasta".format(organism, chain)

modellerpath = "python"#"mod9.19"#"~/.linuxbrew/opt/modeller/bin/mod9.19 "#"~/anaconda3/pkgs/modeller-9.19-py35hddc84ca_1/bin"          #"/usr/local/opt/modeller/bin/"

cmd = "python2 ch_iterate_model.py {0} {1} {2} {3} {4} {5} {6} {7}".format(sequences,   #0
                                                                           workdir,     #1
                                                                           inp_pdb_dir, #2
                                                                           mode,        #3
                                                                           modellerpath,#4
                                                                           num_models,  #5
                                                                           num_proc,    #6
                                                                           chain)
os.system(cmd)
