import csv
import os
import re
import math
from collections import Counter

import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from scipy.stats import entropy


species = ('HomoSapiens', 'MusMusculus')
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


#http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
anchors = {'CDR1': [26, 38], 'CDR2': [55, 65], 'CDR2.5': [80, 86]}
mouse_anchors = {'CDR1': [27, 39], 'CDR2': [56, 66], 'CDR2.5': [81, 88]}


def get_fasta_sequences(inp_file='tcr_aa_sequences.txt', specie='HomoSapiens', chain='A', gene='V', out_file=None):
    if out_file is None:
        out_file = 'TR{}{}_aa_sequences_{}.txt'.format(chain, gene, specie)
    with open(out_file, 'w') as out:
        with open(inp_file) as inp:
            reader = csv.DictReader(inp, delimiter='\t')
            for row in reader:
                if row['#species'] == specie and row['id'].startswith('TR{}{}'.format(chain, gene)):
                    out.write('>{}\n{}\n'.format(row['id'], row['sequence']))


def read_sequences(inp_file, format='fasta'):
    """
    This script will get list of sequences (without ids)
    :param inp_file: str; input file with sequences
    :param format: str; fasta/clustal
    :return: list; list of sequences
    """
    sequences = []
    with open(inp_file, 'r') as inp:
        for record in SeqIO.parse(inp, format):
            sequences.append(record.seq)
    return sequences


def read_sequences_dict(inp_file, format='fasta'):
    """
        This script will get list of sequences
        :param inp_file: str; input file with sequences
        :param format: str; fasta/clustal
        :return: dict; dict of sequences
        """
    sequences = {}
    with open(inp_file, 'r') as inp:
        for record in SeqIO.parse(inp, format):
            sequences[record.id] = str(record.seq)
    return sequences


def get_mafft_alignment(inp_file, outfile, mafft_path='mafft-mac/mafft.bat', params='--reorder --maxiterate 0 --globalpair'):
    if params is None:
        params = ''
    else:
        params += ' '
    os.system('{} {}{} > {}'.format(mafft_path, params, inp_file, outfile))


def get_clustalo_alignment(inp_file, outfile, guidetree=None, force=True):
    cmd = 'clustalo -i {} -o {}'.format(inp_file, outfile)
    if guidetree is not None:
        cmd += ' --guidetree-out {}'.format(guidetree)
    if force is True:
        cmd += ' --force'
    os.system(cmd)


def read_imgt(imgt_file='imgt_all.fasta', specie='HomoSapiens', format='fasta'):
    """
    This script will get dictionary
    :param imgt_file: str; path to file with imgt sequences
    :param specie: str; specie
    :param format: str; format of imgt_file
    :return: dict; dictionary with sequences from imgt: {'gene':sequence}
    """
    imgt_sequences = {}
    specie = specie.lower()
    with open(imgt_file, 'r') as inp:
        for record in SeqIO.parse(inp, format):
            if record.description.split('|')[2].replace(' ', '').lower() == specie:
                imgt_sequences[record.description.split('|')[1]] = record.seq
    return(imgt_sequences)


def get_cdr_posititions(v_gene, specie, imgt=None):
    """
    This script will get dictionary with cdr positions for chosen v_gene.
    :param v_gene: str; identificator of v gene
    :param specie: str; specie
    :param imgt: dict or None; dictionary with tcrs from imgt.
    :return: dict; dictionary with cdr positions: {'CDR1':[start, end],...}
    """
    if imgt is None:
        imgt = read_imgt(specie=specie)
    if specie is 'MusMusculus':
        w_anchors = mouse_anchors
    else:
        w_anchors = anchors
    out_anchors = {}
    imgt_seq = imgt[v_gene]
    for cdr in w_anchors:
        out_anchors[cdr] = []
        imgt_pos = w_anchors[cdr]
        #cdr position in sequence without gaps; +1 and -1 for first gap situation.
        out_anchors[cdr].append(len(str(imgt_seq[:imgt_pos[0]+1]).replace('.', ''))-1)
        out_anchors[cdr].append(len(str(imgt_seq[:imgt_pos[1]]).replace('.', '')))
    return out_anchors


def remove_cdr12(v_gene, tcrs, specie, imgt=None):
    cdr_positions = get_cdr_posititions(v_gene, specie, imgt)
    tcr_sequence = tcrs[v_gene]
    tcr_sequence = tcr_sequence[:cdr_positions['CDR2'][0]] + tcr_sequence[cdr_positions['CDR2'][1]:]
    tcr_sequence = tcr_sequence[:cdr_positions['CDR1'][0]] + tcr_sequence[cdr_positions['CDR1'][1]:]
    return tcr_sequence


def get_pairwise_alignment(tcr_1, tcr_2):
    matrix = matlist.blosum62
    alignment = pairwise2.align.globalds(tcr_1, tcr_2, matrix, min(matrix.values())-1, min(matrix.values())-1)
    return(alignment)
    #print('\n'.join(list('{}\n{}\n{}\n'.format(aln[0], aln[1], aln[2]) for aln in alignment)))


def get_max_score(tcr_1, tcr_2):
    tcr_1_score = get_pairwise_alignment(tcr_1, tcr_1)[0][2]
    tcr_2_score = get_pairwise_alignment(tcr_2, tcr_2)[0][2]
    return max(tcr_1_score, tcr_2_score)


#strange metric, which shows us an angle between two sequences.
#alignment of sequence with itself shows us a length of vector, corresponding to this sequence.
#alignment of two sequences gives us a length of projection of two vectors on surface that is perpendicular to line, that
#connects ends of two sequences...
#added it just in curiosity
def get_max_score_test(tcr_1, tcr_2):
    tcr_1_score = get_pairwise_alignment(tcr_1, tcr_1)[0][2]
    tcr_2_score = get_pairwise_alignment(tcr_2, tcr_2)[0][2]
    aln_score = get_pairwise_alignment(tcr_1, tcr_2)[0][2]
    return math.acos(aln_score/tcr_1_score) + math.acos(aln_score/tcr_2_score)
    #return math.acos(math.sqrt(aln_score**2/(tcr_1_score*tcr_2_score)))


tcrs = read_sequences_dict('position_analysis/TRAV_HomoSapiens.aa_sequences.txt')


def count_column_entropy(inp_file, format='fasta'):
    sequences = read_sequences(inp_file, format)
    lenseq = len(sequences[0])
    num_sequences = len(sequences)
    out_info = []
    for i in range(lenseq):
        column = Counter(seq[i] for seq in sequences)
        probs = list(column[aa] / num_sequences for aa in column)
        entr = entropy(probs)
        gapped = '-' in column
        out_info.append([i, entr, *(column[aa] for aa in amino_acids), gapped])
    out_info = pd.DataFrame(out_info, columns=['#pos', 'entropy', *amino_acids, 'gapped'])
    return out_info


def count_and_write_column_entropy(inp_file, out_file, format):
    out_info = count_column_entropy(inp_file, format)
    out_info.to_csv(out_file, sep='\t', index=None)


def select_positions_align(inpfile, percentile):
    """
    will select best positions with lowest entropy
    :param inpfile: str; output from count_column_entropy
    :param percentile: int; percentile
    :return: list; list of best positions
    """
    inpfile = inpfile[inpfile['gapped']==False]
    inpfile = inpfile.sort_values(by='entropy', axis=0)
    return list(inpfile[:len(inpfile)*percentile//100].sort_values(by='#pos', axis=0)['#pos']) #select percentile


def write_positions_align(inpfile, percentile, outfile):
    positions = select_positions_align(inpfile, percentile)
    with open(outfile, 'w') as out:
        out.write('\t'.join(map(str, positions)))


def select_aminoacid_for_genes(alignment_file, positions, format='fasta'):
    sequences = {}
    prot_aminoacids = {}
    with open(alignment_file, 'r') as inp:
        for record in SeqIO.parse(inp, format):
            sequences[record.id] = str(record.seq)
    for gene in sequences:
        selected_aminoacids = list(sequences[gene][pos] for pos in positions)
        print(selected_aminoacids)
        prot_aminoacids[gene] = selected_aminoacids

    return prot_aminoacids


def select_positions_for_genes(alignment_file, positions, format='fasta'):
    sequences = {}
    prot_positions = {}
    with open(alignment_file, 'r') as inp:
        for record in SeqIO.parse(inp, format):
            sequences[record.id] = str(record.seq)
    for gene in sequences:
        seq_positions = []
        sequence = sequences[gene]
        gaps = list(m.start() for m in re.finditer('-', sequence)) #find positions of all gaps
        gapiti = 0

        for pos in positions:
            while len(gaps) > gapiti and pos > gaps[gapiti]:
                gapiti += 1
            seq_positions.append(pos-gapiti) #change position by gapiti

        prot_positions[gene] = seq_positions

    return prot_positions


def write_positions_for_genes(inpfile, positions, outfile):
    prot_positions = select_positions_for_genes(inpfile,
                               positions)
    with open(outfile, 'w') as out:
        for gene in prot_positions:
            out.write('{}\t{}\n'.format(gene, '\t'.join(map(str, prot_positions[gene]))))


def read_positions_for_genes(inpfile):
    prot_positions = {}
    with open(inpfile, 'r') as inp:
        info = inp.read().split('\n')
        if info[-1] == '':
            info = info[:-1]
        for gene_seq in info:
            prot_positions[gene_seq.split('\t')[0]] = gene_seq.split('\t')[1:]
    return prot_positions


class filename_maker():
    """
    just class for creating fixed names for files
    """
    def __init__(self, chain, gene, specie):
        self.chain = chain
        self.gene = gene
        self.specie = specie
        self.workdir = ''

    def add_dir(self, workdir):
        self.workdir = workdir

    def get_name(self, name):
        return os.path.join(self.workdir, 'TR{}{}_{}.{}.txt'.format(self.chain, self.gene, self.specie, name))


def script_pipeline():
    gene = 'V'
    out_folder = 'position_analysis'
    for spec in species:
        for chain in ['A', 'B']:
            filemaker = filename_maker(chain, gene, spec)
            filemaker.add_dir(out_folder)
            out_file = filemaker.get_name('aa_sequences')
            #os.path.join(out_folder, 'TR{}{}_aa_sequences_{}.txt'.format(chain, gene, spec))
            get_fasta_sequences(specie=spec, chain=chain, gene=gene, out_file=out_file)
            inp_file = out_file
            out_file = filemaker.get_name('align')
            #os.path.join(out_folder, 'TR{}{}_{}.align.txt'.format(chain, gene, spec))
            get_clustalo_alignment(inp_file, out_file, guidetree=filemaker.get_name('guidetree'))
            count_and_write_column_entropy(filemaker.get_name('align'),
                                           filemaker.get_name('entropy'), format='fasta')
            write_positions_for_genes(filemaker.get_name('align'),
                                      select_positions_align(count_column_entropy(filemaker.get_name('align'), 'fasta'),
                                                             20), filemaker.get_name('prot_positions'))

#script_pipeline()
