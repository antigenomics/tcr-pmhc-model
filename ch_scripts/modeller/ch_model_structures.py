from ch_model_base import *

dir_path = os.path.dirname(os.path.realpath(__file__))
home_stdout = sys.stdout

print(sys.argv, len(sys.argv))
if(len(sys.argv) < 5) :
    print("Please provide more arguments")
    sys.exit(0)


template_file = sys.argv[1]
run_mode = sys.argv[2]
num_models = int(sys.argv[3])
chain = sys.argv[4]
if len(sys.argv) > 5:
    inp_pdb_path = sys.argv[5]
else:
    inp_pdb_path = dir_path + '/pdbs/'

structure_name = 'TCR'


template_list = modeller_search_templates(template_file)
print(template_list)
if len(template_list) == 0:
    print('len(template_list) == 0')
    sys.exit(0)

log.verbose()
env = environ()
env.io.atom_files_directory = ['.', inp_pdb_path]
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

aln = alignment(env)

for template_code in template_list:
    mdl = model(env, file=template_code)#, model_segment=('FIRST:', ':LAST'))#('FIRST:', '+{0}:'.format(numaa) + chain))
    aln.append_model(mdl, atom_files=template_code, align_codes=template_code)

aln = modeller_struct_multiple_alignment(aln, 'templates.ali')
aln.append(file='TCR.ali', align_codes='TCR')

aln = modeller_alignment_add_sequence(aln, 'TCR-mult.ali')

f = open('gnuplotfile', 'w')
f.write('plot ')

a = automodel(env, alnfile='TCR-mult.ali',
              knowns=template_list,
              sequence='{0}'.format(structure_name),assess_methods=(assess.DOPE))
a.starting_model = 1
a.ending_model = num_models
a.make()

#profile for template
template_code = template_list[0]
template_mdlch = complete_pdb(env, '{0}'.format(template_code))
alignments_out = 'TCR-mult.ali'

ts = selection(template_mdlch)
ts.assess_dope(output="ENERGY_PROFILE NO_REPORT", file='{0}{1}.profile'.format(structure_name, 0),
               normalize_profile=True, smoothing_window=15)

#profile for structures
for n in range(1, num_models+1):
    mdlch = complete_pdb(env, '{0}.B9999{1:{fill}{align}4}'.format(structure_name, n, fill=0, align='>'))
    s = selection(mdlch)
    s.assess_dope(output="ENERGY_PROFILE NO_REPORT", file='{0}{1}.profile'.format(structure_name, n),
                  normalize_profile=True, smoothing_window=15)

remove_files(num_models)
