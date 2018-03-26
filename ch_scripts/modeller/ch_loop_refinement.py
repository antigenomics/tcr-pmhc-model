from ch_model_base import *

dir_path = os.path.dirname(os.path.realpath(__file__))
home_stdout = sys.stdout

inp_model = sys.argv[1]
chain = sys.argv[2]
num_models = int(sys.argv[3])
cdr3a = sys.argv[4].split('_')
cdr3b = sys.argv[5].split('_')
move_beta = int(sys.argv[6])
inp_pdb_path = sys.argv[7]

if len(sys.argv) == 9:
    structure_name = sys.argv[8]
else:
    structure_name = 'TCR_loop'

log.verbose()
env = environ()
env.io.atom_files_directory = ['.', inp_pdb_path]

class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling

    def select_loop_atoms(self):
        # Two residue ranges (both will be refined simultaneously)
        chain_selection = []
        trim_alpha = 0
        trim_beta = 0
        if chain == 'all' or chain == 'paired' or chain == 'alpha':
            if int(cdr3a[1]) - int(cdr3a[0]) > 6:
                trim_alpha = 3
            print('yay_alpha')
            chain_selection.append(self.residue_range('{0}:A'.format(int(cdr3a[0])+trim_alpha), '{0}:A'.format(int(cdr3a[1])-trim_alpha)))
        if chain == 'all' or chain == 'paired' or chain == 'beta':
            if int(cdr3b[1]) - int(cdr3b[0]) > 6:
                trim_beta = 3
            print('yay_beta')
            chain_selection.append(self.residue_range('{0}:B'.format(int(cdr3b[0])+move_beta+trim_beta), '{0}:B'.format(int(cdr3b[1])+move_beta-trim_beta)))
        print(chain_selection)
        return selection(*chain_selection)

a = MyLoop(env,
           inimodel='{0}'.format(inp_model),  # initial model of the target
           sequence='{0}'.format(structure_name))  # assess each loop with DOPE

a.loop.starting_model = 1           # First loop model
a.loop.ending_model = num_models          # Last loop model
a.loop.md_level = refine.slow

a.make()                            # do modeling and loop refinement

template_mdlch = complete_pdb(env, inp_model)

ts = selection(template_mdlch)
ts.assess_dope(output="ENERGY_PROFILE NO_REPORT", file='{0}{1}.profile'.format(structure_name, 0),
               normalize_profile=True, smoothing_window=15)

for n in range(1, num_models+1):
    mdlch = complete_pdb(env, '{0}.BL{1:{fill}{align}4}0001'.format(structure_name, n, fill=0, align='>'))
    s = selection(mdlch)
    s.assess_dope(output="ENERGY_PROFILE NO_REPORT", file='{0}{1}.profile'.format(structure_name, n),
                  normalize_profile=True, smoothing_window=15)
