import os
import sys
import shutil
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
import analyze_profiles as anpr

BASE_SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))
BASE_HOMEDIR = BASE_SCRIPTDIR[:BASE_SCRIPTDIR.rfind('AdaptiveImm')+len('AdaptiveImm')+1]


def ch_check(info):
    print(info)
    sys.exit()


def crdir(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)


home = os.getcwd()

modeller_home = sys.argv[5]#"/usr/local/opt/modeller/bin/"
if modeller_home == '_':
    modeller_home = ""

sequence_file = sys.argv[1]
workdir = sys.argv[2]
inp_pdb_dir = sys.argv[3]
mode = sys.argv[4]
num_models = sys.argv[6]
num_proc = sys.argv[7]
chain = sys.argv[8]

with open(sequence_file, 'rU') as inp:
    sequences = inp.read().strip().split('>')[1:]

crdir(workdir)
os.chdir(workdir)


def get_struct(struct):
    shutil.copyfile(os.path.join(inp_pdb_dir, 'pdb{0}.ent'.format(struct)), 'pdb{0}.pdb'.format(struct))


def model_struct(num):
    info = sequences[num].split('\n')
    header = info[0] #num|pdbid
    sequence_info = info[1]
    sequence = info[2]

    structures = header.split('|')[1].split(';')
    if chain in ['paired', 'all']:
        move_beta = len(sequence.split('/')[0])
    else:
        move_beta = 0
    pos = header.split('|')[2:-1]
    crdir('seq_'+str(num))
    os.chdir('seq_'+str(num))

    for struct in structures:
        get_struct(struct)

    seq_filename = 'seq.txt'
    with open(seq_filename, 'w') as out:
        out.write('\n'.join(info))

    with open('structure_info.txt', 'w') as out: #name of templates
        out.write('\n'.join(structures))

    with open('TCR.ali', 'w') as out: #alignment file
        out.write(">P1;TCR\nsequence:TCR:::::::0.00: 0.00\n"+str(sequence)+"*\n")

    cmd = "{0} {1}/ch_model_structures.py {2} {3} {4} {5} {6}".format(modeller_home, home,
                                                                          'structure_info.txt',
                                                                          mode, num_models,
                                                                          chain, inp_pdb_dir) #>& model.log #{0}mod9.19

    print(cmd)
    #os.system(cmd)
    profile = int(anpr.get_profile_array('', seq_filename, chain, int(num_models), 'TCR'))
    #inp_model = 'TCR.B9999{0:{fill}{align}4}'.format(profile, fill=0, align='>')
    inp_model = 'pdb1oga.pdb'
    num_loop_model = 3
    if chain in ['paired', 'all']:
        cmd = "{0} {1}/ch_loop_refinement.py {2} {3} {4} {5} {6} {7} {8} {9}".format(modeller_home, home, inp_model, 'alpha', num_loop_model, pos[4], pos[5], move_beta, '.', 'TCR_loop_alpha')
        os.system(cmd)
        profile = anpr.get_profile_array('', seq_filename, chain, num_loop_model, 'TCR_loop_alpha', use_default_struct=False, cdr3only=True)
        inp_model = '{0}.BL{1:{fill}{align}4}0001'.format('TCR_loop_alpha', profile, fill=0, align='>')
        cmd = "{0} {1}/ch_loop_refinement.py {2} {3} {4} {5} {6} {7} {8} {9}".format(modeller_home, home, inp_model, 'beta', num_loop_model, pos[4], pos[5], move_beta, '.', 'TCR_loop_beta')
        os.system(cmd)
        profile = anpr.get_profile_array('', seq_filename, chain, num_loop_model, 'TCR_loop_beta', use_default_struct=False,
                                         cdr3only=True)
        print(profile)
        profile = anpr.get_profile_array('', seq_filename, chain, num_loop_model, 'TCR_loop_beta', use_default_struct=True,
                                         cdr3only=True)
        print(profile)
    else:
        cmd = "{0} {1}/ch_loop_refinement.py {2} {3} {4} {5} {6} {7} {8} {9}".format(modeller_home, home, inp_model,
                                                                                 chain, num_loop_model, pos[4], pos[5], move_beta,
                                                                                 '.', 'TCR_loop')
        os.system(cmd)
        profile = anpr.get_profile_array('', seq_filename, chain, num_loop_model, 'TCR_loop', use_default_struct=False,
                                         cdr3only=True)
        print(profile)
        profile = anpr.get_profile_array('', seq_filename, chain, num_loop_model, 'TCR_loop', use_default_struct=True,
                                         cdr3only=True)
        print(profile)
    os.chdir('..')

    '''
    numfiles = (len(os.listdir('.')))
    
    if numfiles < 6:
        os.chdir('..')
        shutil.rmtree('seq_'+str(num))
    else:
        profile = int(anpr.get_profile_array('', type))
        if profile == 10:
            outfile = 'TCR.B99990010.pdb'
        else:
            outfile = 'TCR.B9999000{0}.pdb'.format(profile)
        shutil.copyfile(outfile, '../model_structs/seq_{0}.pdb'.format(num))
        for fitmodel in glob.glob('*fit*'):
            os.remove(fitmodel)
        os.chdir('..')
    '''

crdir('model_structs')
model_struct(0)
#pool = Pool(num_proc)
#pool.map(model_struct, [k for k in range(len(sequences))])
#pool.close()
#pool.join()