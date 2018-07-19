import csv
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd


#those nucleotides will be added to J segment when function translate_aa_seq is used.
two_c_nucleotides_dict = {"homosapiens":{"TRAC": "AT", "TRBC": "AG", "TRGC": "AT", "TRDC": "GA"},
                     "musmusculus":{"TRAC": "AC", "TRBC": "AG", "TRGC": "AC", "TRDC": "AA"},
                     "macacamulatta":{"TRAC" :"AT", "TRBC": "AG", "TRGC": "!!!!!!!!!", "TRDC": "AA"}} #will cause errors. There is no TRGC for macacamulatta in imgt


#get aminoacid sequence for V/J
def translate_aa_seq(nucleotide_sequence, reference_point, segment_type, specie, gene):
    try:
        if segment_type in ('J', 'Joining'):
            coding_dna = Seq(nucleotide_sequence[reference_point+1:] +
                             two_c_nucleotides_dict[specie.lower()]["{}C".format(gene[:3])],
                             generic_dna)
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
                row['sequence'] = translate_aa_seq(row['sequence'], int(row['reference_point']), row['segment'],
                                                   row['#species'], row["id"])

                for pos in cdr_names:
                    if int(row[pos]) > -1:
                        row[pos] = int(row[pos])//3
                #genes from MacacaMulatta sometimes have problems here: TRBV7-2*01
                writer.writerow({key: row[key] for key in out_columns})


class tcr_constructor():
    def __init__(self,
                 vj_sequences_file,
                 c_sequences_file,
                 segments_file='segments.aaparts.txt'):
        self.tcr_dict = None
        self.c_dict = None
        self.vj_sequences_file = vj_sequences_file
        self.c_sequences_file = c_sequences_file
        self.segments = None
        self.segments_file = segments_file

    def get_tcr_dict(self):
        return self.tcr_dict

    def get_c_dict(self):
        return self.c_dict

    # set tcr_dict
    def set_tcr_genes(self):
        if self.tcr_dict is not None:
            return

        self.tcr_dict = {}

        with open(self.vj_sequences_file, 'r') as inp:
            reader = csv.DictReader(inp, delimiter='\t')
            for row in reader:
                if row['#species'] not in self.tcr_dict:
                    self.tcr_dict[row['#species']] = {}
                self.tcr_dict[row['#species']][row['id']] = row['sequence']

    #set c_dict
    def set_c_genes(self):
        if self.c_dict is not None:
            return
        self.c_dict = {}
        with open(self.c_sequences_file, 'r') as inp:
            reader = csv.DictReader(inp, delimiter='\t')
            for row in reader:
                if row['species'] not in self.c_dict:
                    self.c_dict[row['species']] = {}
                self.c_dict[row['species']][row['c.id']] = row['c.seq']

    def get_segments(self):
        return self.segments

    # set segments from https://github.com/antigenomics/vdjdb-db/blob/master/res/segments.aaparts.txt
    def set_segment_prediction(self):
        if self.segments is not None:
            return
        self.segments = pd.read_table(self.segments_file)

    # predict v/j genes using cdr3 sequence
    def predict_segment(self, cdr3, chain_type, gene_type, specie):
        if chain_type is None:
            raise TypeError('chain_type must be str')
        self.set_segment_prediction()
        w_segments = self.segments[(self.segments['species'] == specie) & (self.segments['gene'] == chain_type) & (
                    self.segments['type'] == gene_type)]
        while len(cdr3) > 0:
            if len(w_segments[w_segments['cdr3'] == cdr3]) > 0:
                return str(w_segments[w_segments['cdr3'] == cdr3]['segm'].reset_index(drop=True)[0]) + '*01'
            if gene_type == 'J':
                cdr3 = cdr3[1:]
            elif gene_type == 'V':
                cdr3 = cdr3[:-1]
        return None


class tcr_aminoacid_constructor(tcr_constructor):
    def __init__(self,
                 vj_sequences_file='tcr_aa_sequences.txt',
                 c_sequences_file='constant.txt',
                 segments_file='segments.aaparts.txt'):
        tcr_constructor.__init__(self, vj_sequences_file, c_sequences_file)


    def get_tcr(self, cdr3, specie, v_gene=None, j_gene=None, add_c=False,
                try_to_predict_missing=False, chain_type=None, force=False, seqlengths=False):
        """
        get tcr sequence based on specie, cdr3, v_gene and j_genes
        :param cdr3: str; cdr3 sequence
        :param specie: str; specie name (HomoSapiens, MacacaMulatta, MusMusculus)
        :param v_gene: str; V segment with *01 at the end
        :param j_gene: str; J segment with *01 at the end
        :param add_c: True/False; if True, constant region will be added
        :param try_to_predict_missing: True/False; if True, then missing v_gene and j_gene
        will be predicted based on cdr3 sequence. if True, then chain_type is needed.
        :param chain_type: str; TRA/TRB
        :param force: TRUE/FALSE; if True, then tcrs without v_gene won't result in stopped calculations; tcrs with missing
        v segment will have "NO_V_SEGMENT_FOUND" sequences
        :param seqlengths: TRUE/FALSE; if True, then in return you will get tuple with tcr sequence and lengths of v_gene, cdr3, j_gene and c_gene
        :return: str; tcr sequence
        """
        self.set_tcr_genes()

        if v_gene is None:
            if try_to_predict_missing is True:
                v_gene = self.predict_segment(cdr3, chain_type, 'V', specie)
            if v_gene is None and force is True:
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
                c_seq = self.c_dict[specie]['TRBC{}*01'.format(j_gene[4])]
            else:
                c_seq = self.c_dict[specie]['TRAC*01']

        ret_tcr = (self.tcr_dict[specie][v_gene][:-1], cdr3, j_seq, c_seq)
        if seqlengths is False:
            return "".join(ret_tcr)
        elif seqlengths is True:
            return "".join(ret_tcr), (*map(len, ret_tcr))

#USAGE:
#tcr_maker = tcr_constuctor()
#
#print(tcr_maker.get_tcr(cdr3='CVVNSPNDYKLSF', specie='HomoSapiens', add_c=True, try_to_predict_missing=True, chain_type='TRA', force=True))

class tcr_nucleotide_constructor(tcr_constructor):
    def __init__(self,
                 vj_sequences_file='tcr_genes.txt',
                 c_sequences_file='constant.nucleotide.txt'):
        tcr_constructor.__init__(self, vj_sequences_file, c_sequences_file)
        self.cdr3points = None

    def set_tcr_genes(self):
        if self.tcr_dict is not None:
            return

        self.tcr_dict = {}

        with open(self.vj_sequences_file, 'r') as inp:
            reader = csv.DictReader(inp, delimiter='\t')
            for row in reader:
                if row['#species'] not in self.tcr_dict:
                    self.tcr_dict[row['#species']] = {}
                self.tcr_dict[row['#species']][row['id']] = (row['sequence'], int(row['reference_point']))

    def get_cdr3_nucleotide_sequence(self, cdr3):
        return "N"*len(cdr3)*3

    def get_tcr(self, cdr3, specie, cdr3vend=1, cdr3jstart=0, v_gene=None, j_gene=None, add_c=False,
                try_to_predict_missing=False, chain_type=None, force=False, seqlengths=False):
        """
        get tcr sequence based on specie, cdr3, v_gene and j_genes
        :param cdr3: str; cdr3 aminoacid sequence
        :param specie: str; specie name (HomoSapiens, MacacaMulatta, MusMusculus)
        :param cdr3vend: int; part of v-gene, that is used in cdr3 plus one (in aminoacids)
        :param cdr3jstart: int; part of j-gene, that is used in cdr3 (in aminoacids)
        :param v_gene: str; V segment with *01 at the end
        :param j_gene: str; J segment with *01 at the end
        :param add_c: True/False; if True, constant region will be added
        :param try_to_predict_missing: True/False; if True, then missing v_gene and j_gene
        will be predicted based on cdr3 sequence. if True, then chain_type is needed.
        :param chain_type: str; TRA/TRB
        :param force: TRUE/FALSE; if True, then tcrs without v_gene won't result in stopped calculations; tcrs with missing
        v segment will have "NO_V_SEGMENT_FOUND" sequences
        :param seqlengths: TRUE/FALSE; if True, then in return you will get tuple with tcr sequence and lengths of v_gene, cdr3, j_gene and c_gene
        :return: str; tcr sequence
        """
        self.set_tcr_genes()


        if cdr3vend < 0:
            cdr3vend = 1
        if v_gene is None:
            if try_to_predict_missing is True:
                v_gene = self.predict_segment(cdr3, chain_type, 'V', specie)
            if v_gene is None and force is True:
                return "NO_V_SEGMENT_FOUND"
        v_info = self.tcr_dict[specie][v_gene]
        v_point = v_info[1]
        v_seq = v_info[0][:v_point+3*(cdr3vend-1)]

        j_seq = ''
        if j_gene is None and try_to_predict_missing is True:
            j_gene = self.predict_segment(cdr3, chain_type, 'J', specie)

        if j_gene is not None:
            j_info = self.tcr_dict[specie][j_gene]
            j_point = j_info[1]
            j_seq = j_info[0][(j_point+1)-3*(len(cdr3) - cdr3jstart - 1):] #(-1) - for last F;
            #print(j_info[0], j_point, j_point+1, 3*(len(cdr3)- cdr3jstart-1))

        c_seq = ''
        if add_c is True:
            self.set_c_genes()
            if j_gene[2] == 'B':
                c_seq = self.c_dict[specie]['TRBC{}*01'.format(j_gene[4])]
            else:
                c_seq = self.c_dict[specie]['TRAC*01']

        ret_tcr = (v_seq, self.get_cdr3_nucleotide_sequence(cdr3[cdr3vend:cdr3jstart]), j_seq, c_seq[1:])
        #print(ret_tcr)
        if seqlengths is False:
            return "".join(ret_tcr)
        elif seqlengths is True:
            return "".join(ret_tcr), (*map(len, ret_tcr))


#tcr_n_constructor = tcr_nucleotide_constructor()
#tcr_a_constructor = tcr_aminoacid_constructor()
#tcr_construct = tcr_n_constructor.get_tcr("CASSLAPGATNEKLFF", "HomoSapiens", 5, 9, "TRBV7-6*01", "TRBJ1-4*01", add_c=True, seqlengths=True)
#tcr_aminoacid = tcr_a_constructor.get_tcr("CASSLAPGATNEKLFF", "HomoSapiens", "TRBV7-6*01", "TRBJ1-4*01", add_c=True, seqlengths=True)
#
#print(tcr_construct[0])
#print(Seq(tcr_construct[0], generic_dna).translate())
#print(tcr_aminoacid[0])