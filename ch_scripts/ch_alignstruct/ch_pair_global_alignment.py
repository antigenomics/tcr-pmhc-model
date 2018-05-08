import os
import pandas as pd
import numpy as np
import argparse

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will get information from file with distances between "
                                             "superimposed PDB structures and then will get pairwise alignments. \
Current version is %s" % curr_version)

parser.add_argument("-i", help="Path to file with distances between structures. Default is "
                               "'../seqdata/HomoSapiens/distances/distances_from_pdb_BCDR3.txt'", required=False, dest="input_dir")
parser.add_argument("-o", help="Path to output directory. "
                               "Default is '../seqdata/HomoSapiens/distances/pairwise_alignments_BCDR3.txt'",
                    required=False, dest='output_dir')

myargs = parser.parse_args()

if myargs.input_dir is None:
    myargs.input_dir = 'seqdata/HomoSapiens/distances/distances_from_pdb_BCDR2.txt'
if myargs.output_dir is None:
    myargs.output_dir = 'seqdata/HomoSapiens/distances/pairwise_alignments_BCDR2.txt'

def traceback_prc(wmatrix, iseq, jseq, gapmatrix):

    def funcname(gdict, ch1, ch2, gpscore):
        nch = gdict['way']
        gscore = 0
        mindistindex = gdict['distances'].index(min(gdict['distances']))
        for k in range(len(gdict['distances'])):
            kk = list(map(int, gdict['gaps'][k].split(';')))
            if k == mindistindex:
                ch1 += iseq[kk[0]]
                ch2 += jseq[kk[1]]
            else:
                if nch == 'i':
                    ch1 += '-'
                    ch2 += jseq[kk[1]]
                    gscore += gdict['distances'][k]
                else:
                    ch1 += iseq[kk[0]]
                    ch2 += '-'
                    gscore += gdict['distances'][k]
        gpscore[0] += gscore
        gpscore[1] += len(gdict['distances'])-1
        return ch1, ch2, gpscore

    i = len(iseq) - 1
    j = len(jseq) - 1
    ioutseq = ''
    joutseq = ''
    father = wmatrix[i][j]['father'].split(';')
    gapdict = {'way':'None', 'gaps':[], 'distances':[]}
    gapscore = [0, 0]
    while father != ['None']:
        father = list(map(int, father))
        if father[0] == int(i)-1 and father[1] == int(j)-1:
            if gapdict['way'] != 'None':
                gapdict['gaps'].append(str(i) + ';' + str(j))
                gapdict['distances'].append(gapmatrix[j][i])
                ioutseq, joutseq, gapscore = funcname(gapdict, ioutseq, joutseq, gapscore)
                gapdict = {'way': 'None', 'gaps': [], 'distances': []}
            else:
                ioutseq += iseq[i]
                joutseq += jseq[j]
        elif father[0] == int(i):
            if gapdict['way'] == 'j':
                ioutseq, joutseq, gapscore = funcname(gapdict, ioutseq, joutseq, gapscore)
            else:
                gapdict['way'] = 'i'
                gapdict['gaps'].append(str(i)+';'+str(j))
                gapdict['distances'].append(gapmatrix[j][i])
        else:
            if gapdict['way'] == 'i':
                ioutseq, joutseq, gapscore = funcname(gapdict, ioutseq, joutseq, gapscore)
            else:
                gapdict['way'] = 'j'
                gapdict['gaps'].append(str(i) + ';' + str(j))
                gapdict['distances'].append(gapmatrix[j][i])
        i = father[0]
        j = father[1]
        father = wmatrix[i][j]['father'].split(';')
    if gapdict['way'] != 'None':
        gapdict['gaps'].append(str(i) + ';' + str(j))
        gapdict['distances'].append(gapmatrix[j][i])
        ioutseq, joutseq, gapscore = funcname(gapdict, ioutseq, joutseq, gapscore)
        gapdict = {'way': 'None', 'gaps': [], 'distances': []}
    else:
        ioutseq += iseq[i]
        joutseq += jseq[j]
    if gapscore[1] != 0:
        gapscore = (gapscore[0]/gapscore[1])
    else:
        gapscore = 0
    return (str(ioutseq)[::-1]), (str(joutseq)[::-1]), gapscore

stralignment = pd.read_csv(myargs.input_dir, sep='\t')

with open(myargs.output_dir, 'w') as out:
    out.writelines('pdb1\tpdb2\tscore\tgapscore\tseq1\tseq2\n')
    pdbs = stralignment.groupby("pdb1")
    for pdb1 in pdbs:
        for pdb2 in pdb1[1].groupby("pdb2"):
            pdmatrix = pdb2[1].pivot('#1', '#2', 'distance')
            print()
            seq1 = pdb2[1]['pdb1seq'][pdb2[1]['pdb1seq'].index[0]]
            seq2 = pdb2[1]['pdb2seq'][pdb2[1]['pdb2seq'].index[0]]
            pdmatrix.index.name = None
            pdmatrix.columns.name = None
            shapes = (np.shape(pdmatrix))
            workmatrix = [[0 for j in range(shapes[1])] for i in range(shapes[0])]
            i = j = 0
            while i < shapes[0] and j < shapes[1]:
                pdi = j
                pdj = i
                stepi = stepj = 1
                listmatrices = []
                if i == 0:
                    stepi = 0
                else:
                    listmatrices.append(workmatrix[i - stepi][j])
                if j == 0:
                    stepj = 0
                else:
                    listmatrices.append(workmatrix[i][j - stepj])
                if i == 0 and j == 0:
                    workmatrix[i][j] = {'distance': -pdmatrix[pdi][pdj], 'father': 'None', 'ij': str(i) + ';' + str(j)}
                else:
                    listmatrices.append(workmatrix[i - stepi][j - stepj])
                    values = [listmatrices[k]['distance'] for k in range(len(listmatrices))]
                    minmatrix = listmatrices[values.index(max(values))]
                    workmatrix[i][j] = {'distance': minmatrix['distance'] - pdmatrix[pdi][pdj], 'father':minmatrix['ij'], 'ij':str(i)+';'+str(j)}
                i += 1
                if i == shapes[0]:
                    i = 0
                    j += 1
            print(pdb1[0], pdb2[0])
            score = (-workmatrix[-1][-1]['distance']/max(len(seq1), len(seq2)))
            alseq1, alseq2, gapscore = traceback_prc(workmatrix, seq1, seq2, pdmatrix)
            out.writelines('\t'.join([str(pdb1[0]), str(pdb2[0]), str(score), str(gapscore), str(alseq1), str(alseq2)])+'\n')