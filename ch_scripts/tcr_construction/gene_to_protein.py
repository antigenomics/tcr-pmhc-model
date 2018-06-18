import csv
import warnings

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd


#get aminoacid sequence for V/J
def translate_aa_seq(nucleotide_sequence, reference_point, segment_type):
    try:
        if segment_type in ('J', 'Joining'):
            coding_dna = Seq(nucleotide_sequence[reference_point+1:], generic_dna)
        elif segment_type in ('V', 'Variable'):
            coding_dna = Seq(nucleotide_sequence[:reference_point], generic_dna)
        else:
            raise UserWarning
    except UserWarning:
        warnings.warn('your segment_type ({}) is not V/Variable nor J/Joining. Full aminoacid sequence will be returned'.format(segment_type))
        coding_dna = Seq(nucleotide_sequence, generic_dna)
    return str(coding_dna.translate())


#read file with nucleotide sequences and reference points and then write file with aminoacid sequences
def write_aa_sequences(inpfile='tcr_genes.txt', outfile='tcr_aa_sequences.txt', use_cdr_names=True):
    out_columns = ['#species', 'gene', 'segment', 'id', 'sequence']
    if use_cdr_names is True:
        cdr_names = ['cdr{}.{}'.format(i, k) for i in (1, 2, 2.5) for k in ('start', 'end')]
        out_columns += list(cdr_names)
    else:
        cdr_names = []

    with open(outfile, 'w') as out:
        writer = csv.DictWriter(out, fieldnames=out_columns, delimiter='\t')
        writer.writeheader()
        with open(inpfile, 'r') as inp:
            reader = csv.DictReader(inp, delimiter='\t')
            for row in reader:
                if row['#species'] not in ('MusMusculus', 'HomoSapiens'):
                    continue
                row['sequence'] = translate_aa_seq(row['sequence'], int(row['reference_point']), row['segment'])

                for pos in cdr_names:
                    if int(row[pos]) > -1:
                        row[pos] = int(row[pos])//3
                #genes from MacacaMulatta sometimes have problems here: TRBV7-2*01
                writer.writerow({key: row[key] for key in out_columns})


class tcr_constuctor():
    def __init__(self):
        self.tcr_dict = None
        self.c_dict = None
        self.segments = None

    def get_tcr_dict(self):
        return self.tcr_dict

    def get_c_dict(self):
        return self.c_dict

    def get_segments(self):
        return self.segments

    #set tcr_dict
    def set_tcr_genes(self, inpfile='tcr_aa_sequences.txt'):
        if self.tcr_dict is not None:
            return

        self.tcr_dict = {}

        with open(inpfile, 'r') as inp:
            reader = csv.DictReader(inp, delimiter='\t')
            for row in reader:
                if row['#species'] not in self.tcr_dict:
                    self.tcr_dict[row['#species']] = {}
                self.tcr_dict[row['#species']][row['id']] = row['sequence']

    #set c_dict
    def set_c_genes(self, inpfile='constant.txt'):
        if self.c_dict is not None:
            return
        self.c_dict = {}
        with open(inpfile, 'r') as inp:
            reader = csv.DictReader(inp, delimiter='\t')
            for row in reader:
                if row['species'] not in self.c_dict:
                    self.c_dict[row['species']] = {}
                self.c_dict[row['species']][row['c.id']] = row['c.seq']

    #set segments from https://github.com/antigenomics/vdjdb-db/blob/master/res/segments.aaparts.txt
    def set_segment_prediction(self, inpfile='segments.aaparts.txt'):
        if self.segments is not None:
            return
        self.segments = pd.read_table(inpfile)

    #predict v/j genes using cdr3 sequence
    def predict_segment(self, cdr3, chain_type, gene_type, specie):
        if chain_type is None:
            raise TypeError('chain_type must be str')
        self.set_segment_prediction()
        w_segments = self.segments[(self.segments['species'] == specie) & (self.segments['gene'] == chain_type) & (self.segments['type'] == gene_type)]
        while len(cdr3) > 0:
            if len(w_segments[w_segments['cdr3'] == cdr3]) > 0:
                return str(w_segments[w_segments['cdr3'] == cdr3]['segm'].reset_index(drop=True)[0]) + '*01'
            if gene_type == 'J':
                cdr3 = cdr3[1:]
            elif gene_type == 'V':
                cdr3 = cdr3[:-1]
        return None

    def get_tcr(self, cdr3, specie, v_gene=None, j_gene=None, add_c=False,
                try_to_predict_missing=False, chain_type=None, force=False):
        """
        get tcr sequence based on specie, cdr3, v_gene and j_genes
        :param cdr3: str; cdr3 sequence
        :param specie: str; specie name (HomoSapiens, MacacaMulatta, MusMusculus)
        :param v_gene: str; V segment
        :param j_gene: str; J segment
        :param add_c: True/False; if True, constant region will be added
        :param try_to_predict_missing: True/False; if True, then missing v_gene and j_gene
        will be predicted based on cdr3 sequence. if True, then chain_type is needed.
        :param chain_type: str; TRA/TRB
        :param force: TRUE/FALSE; if True, then tcrs without v_gene won't stop calculations; tcrs with missing
        v segment will have "NO_V_SEGMENT_FOUND" sequences
        :return: str; tcr sequence
        """
        self.set_tcr_genes()

        if v_gene is None:
            if try_to_predict_missing is True:
                v_gene = self.predict_segment(cdr3, chain_type, 'V', specie)
            if v_gene is None and force == True:
                return "NO_V_SEGMENT_FOUND"

        j_seq = ''
        if j_gene is None and try_to_predict_missing is True:
            j_gene = self.predict_segment(cdr3, chain_type, 'J', specie)

        if j_gene is not None:
            j_seq = self.tcr_dict[specie][j_gene][1:]

        c_seq = ''
        if add_c is True:
            self.set_c_genes()
            if j_gene[2] == 'B':
                c_seq = self.c_dict[specie]['TRBC{}'.format(j_gene[4])]
            else:
                c_seq = self.c_dict[specie]['TRAC']

        return self.tcr_dict[specie][v_gene][:-1] + cdr3 + j_seq + c_seq


tcr_maker = tcr_constuctor()
print(tcr_maker.get_tcr(cdr3='CVVNSPNDYKLSF', specie='HomoSapiens', add_c=True, try_to_predict_missing=True, chain_type='TRA', force=True))
#def get_tcr_genes(inpfile='tcr_aa_sequences.txt'):
#    tcr_dict = {}
#    with open(inpfile, 'r') as inp:
#        reader = csv.DictReader(inp, delimiter='\t')
#        for row in reader:
#            if row['#species'] not in tcr_dict:
#                tcr_dict[row['#species']] = {}
#            tcr_dict[row['#species']][row['id']] = row['sequence']
#    return(tcr_dict)
#
#
#tcr_dict = get_tcr_genes()


#def get_tcr(v_gene, j_gene, cdr3, species, tcr_dict=tcr_dict, add_c=False):
#    if j_gene == '':
#        return(tcr_dict[species][v_gene][:-1]+cdr3)
#    c_gene = ''
#    if add_c is True:
#
#    else:
#        return(tcr_dict[species][v_gene][:-1]+cdr3+tcr_dict[species][j_gene][1:]+c_gene)

