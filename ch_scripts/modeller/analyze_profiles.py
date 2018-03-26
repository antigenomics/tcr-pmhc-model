import os
import pandas as pd
from Bio.PDB.Polypeptide import three_to_one


def get_id_from_fasta(filename):
    inp = open(filename, 'r')
    info = inp.readline().replace('>', '', 1).strip().split('|')
    inp.close()
    if len(info) == 8:
        return info[1]
    else:
        return info[-2]


def get_fasta_info(filename, type):
    """
    This function will get info from header of FASTA sequence from file (now unused)
    :param filename: input file in FASTA format
    :param type: type of file ('seq' or 'pdb')
    :return: list with information from file
    """
    inp = open(filename, 'r')
    info = inp.readline().replace('>', '', 1).strip().split('|')
    inp.close()
    if type == 'pdb':
        return(['pdb', info[3], info[4], info[5]])
    else:
        return(['justseq', info[5], info[6], info[7]])


def get_cdrs(seqinfo='seq', inppath='', type='seq'):
    """
    This function is used to get CDR sequences from file in FASTA format for blast (now unused)
    :param seqinfo: file in FASTA format
    :param inppath: folder with file in FASTA format
    :param type: type of file ('seq' or 'pdb'). On it depends what info you can get from header
    :return: dictionary with keys: 'CDR1', 'CDR2', 'CDR3'
    """
    info = get_fasta_info('{0}'.format(os.path.join(inppath, seqinfo)), type)
    cdrs = {'CDR1': info[1], 'CDR2': info[2], 'CDR3': info[3]}
    return(cdrs)


def get_profile(inppath, profilename, seq=''):
    """
    This function will read profile file with DOPE energies and return DataFrame based on information from it
    :param inppath: path to profile (with profile included)
    :param profilename: name of the profile (needed for DataFrame)
    :param seq: just sequence on which structure for profile is based (if '', then it will be computed from profile)
    :return: DataFrame and seq
    """
    energy = []
    pos = []
    aa = []
    with open(inppath, 'r') as inp:
        for line in inp:
            if not line.startswith('#') and len(line) > 10:
                info = line.strip().split()

                pos.append(int(info[0].strip()))
                energy.append(float(info[-1]))
                if info[1] == 'UNK':
                    aa.append('X')
                else:
                    aa.append(three_to_one(info[1]))
        energy_profile = pd.DataFrame.from_items(
            [('pos', pos), ('aa', aa), ('profile', [str(profilename)] * len(pos)), ('energy', energy)])

        if seq == '':
            seq = ''.join(aa)
        else:
            assert seq == ''.join(seq)
    return energy_profile, seq


def get_tcr_info(inppath, chains):
    """
    This function will get you info from file in PIR-like format and return positions for CDRs
    :param inppath: input file in PIR-like format
    :param chains: alpha/beta/paired(=all)
    :return: CDR positions
    """
    with open(inppath, 'r') as inp:
        info = inp.read().split('\n')
        sequences = info[2].split('/')
        cdr_positions = [map(int, pos.split('_')) for pos in info[0].split('|')[2:-1]]
    if chains in ['paired', 'all']:
        return [cdr_positions[i] for i in (0, 2, 4)] + \
               [[cdr_positions[i][0]+len(sequences[0]), cdr_positions[i][1]+len(sequences[0])] for i in (1, 3, 5)]
    if chains == 'alpha':
        return [cdr_positions[i] for i in (0, 2, 4)]
    elif chains == 'beta':
        return [cdr_positions[i] for i in (1, 3, 5)]

def get_profile_array(inppath, filename, chains, num=10, profilename = 'TCR', use_default_struct=False, cdr3only=False):
    cdrs = get_tcr_info(os.path.join(inppath, filename), chains)
    energy_profiles = pd.DataFrame(columns=['pos', 'aa', 'profile', 'energy'])
    seq = ''
    if use_default_struct is True:
        w_n = 1
        w_k = 0
    else:
        w_n = 0
        w_k = 1
    for i in range(num+w_n):
        profile = i+w_k
        energy_profile, seq = get_profile('{0}{1}{2}.profile'.format(inppath, profilename, i+w_k), profile, seq)
        energy_profiles = pd.concat([energy_profiles, energy_profile], axis=0)

    energy_profiles['mean'] = energy_profiles['energy'].groupby(energy_profiles['pos']).transform(lambda x: x.mean())
    energy_profiles['std'] = energy_profiles['energy'].groupby(energy_profiles['pos']).transform(lambda x: x.std())
    if len(cdrs) == 6:
        beta_trim = 0 #trim cdr3beta
        alpha_trim = 0 #trim cdr3alpha
        if cdrs[5][1]-cdrs[5][0] > 6:
            beta_trim = 3
        if cdrs[4][1]-cdrs[4][0] > 6:
            alpha_trim = 3

        if cdr3only:
            energy_profiles['sumECDRs'] = energy_profiles['energy'].groupby(
                energy_profiles['profile']).transform(
                lambda x: x.loc[cdrs[4][0] + alpha_trim:cdrs[4][1] - alpha_trim].sum() / (
                        cdrs[4][1] - cdrs[4][0] - 2 * alpha_trim) +
                          x.loc[cdrs[5][0] + beta_trim:cdrs[5][1] - beta_trim].sum() / (
                                  cdrs[5][1] - cdrs[5][0] - 2 * beta_trim))
        else:
            energy_profiles['sumECDRs'] = energy_profiles['energy'].groupby(
                energy_profiles['profile']).transform(
                lambda x: x.loc[cdrs[4][0] + alpha_trim:cdrs[4][1] - alpha_trim].sum() / (
                        cdrs[4][1] - cdrs[4][0] - 2 * alpha_trim) +
                          x.loc[cdrs[5][0] + beta_trim:cdrs[5][1] - beta_trim].sum() / (
                                  cdrs[5][1] - cdrs[5][0] - 2 * beta_trim) +
                          sum([x.loc[cdrs[i][0]:cdrs[i][1]].sum()/(cdrs[i][1]-cdrs[i][0]) for i in (0, 1, 2, 3)])) #cdr1 a/b, cdr2 a/b
    else:
        trim = 0
        if cdrs[2][1]-cdrs[2][0] > 6:
            trim = 3
        if cdr3only:
            energy_profiles['sumECDRs'] = energy_profiles['energy'].groupby(
                energy_profiles['profile']).transform(
                lambda x: x.loc[cdrs[2][0] + trim:cdrs[2][1] - trim].sum() / (cdrs[2][1] - cdrs[2][0] - 2 * trim))
        else:
            energy_profiles['sumECDRs'] = energy_profiles['energy'].groupby(
                energy_profiles['profile']).transform(
                lambda x: x.loc[cdrs[2][0] + trim:cdrs[2][1] - trim].sum() / (cdrs[2][1] - cdrs[2][0] - 2 * trim) +
                          sum([x.loc[cdrs[i][0]:cdrs[i][1]].sum() / (cdrs[i][1] - cdrs[i][0]) for i in (0, 1)])) #cdr1, cdr2

    profile = (energy_profiles.groupby(['sumECDRs']).min().iloc[0]['profile'])
    energy_profiles.to_csv('{0}energy_profiles.txt'.format(inppath), sep='\t', index=None)
    return profile
    #energy_profiles_unmelted = energy_profiles.pivot_table(index=['pos', 'aa', 'std', 'mean'], columns='profile', values='energy').reset_index().sort_values(by=['pos'])
    #energy_profiles_unmelted.columns.name = None
    #print(energy_profiles)
    #print(energy_profiles_unmelted.loc[cdrpos['CDR3'][0]+3:cdrpos['CDR3'][1]-3]['mean'].sum()/(len(cdrpos['CDR3'])-6))

#for i in range(18):
#    inppath = '{0}/seq_{1}/'.format('musmusculus_beta_sequences_test_2', i)
#    profile = int(get_profile_array(inppath))
#    print(profile)

#Code for some tests..
#pdbname = get_id_from_fasta('{0}seq'.format(inppath))
#real_profile, _ = get_profile('{0}{1}.profile'.format(inppath, pdbname), pdbname)
#
#print(real_profile.loc[cdrpos['CDR3'][0]+3:cdrpos['CDR3'][1]-3]['energy'].sum()/(len(cdrpos['CDR3'])-6) + real_profile.loc[cdrpos['CDR2'][0]:cdrpos['CDR2'][1]]['energy'].sum()/(len(cdrpos['CDR2'])) + real_profile.loc[cdrpos['CDR1'][0]:cdrpos['CDR1'][1]]['energy'].sum()/(len(cdrpos['CDR1'])))
#energy_profiles = energy_profiles.merge(real_profile, left_on=['pos', 'aa'], right_on=['pos', 'aa'])
#energy_profiles['lowertsigma'] = (energy_profiles['mean'] - energy_profiles['std']*3 < energy_profiles['energy'])
#energy_profiles['uppertsigma'] = (energy_profiles['mean'] + energy_profiles['std']*3 > energy_profiles['energy'])
#print(energy_profiles.loc[cdrpos['CDR3'][0]:cdrpos['CDR3'][1]])
