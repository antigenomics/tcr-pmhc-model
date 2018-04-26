from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import csv
import warnings

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


#def get_cdrs(inpfile='aa_cdrs.txt'):
#    cdrs = {}
#    with open(inpfile, 'r') as inp:
#        reader = csv.DictReader(inp, delimiter='\t')
#        for row in reader:
#            if row['species'] not in cdrs:
#                cdrs[row['species']] = {}
#            cdrs[row['species']][row['gene']] = {}
#            for pos in ['cdr1aa', 'cdr2aa', 'cdr2.5aa']:
#                cdrs[row['species']][row['gene']][pos]=row[pos]
#    return cdrs
#
#
#test_cdrs =  get_cdrs()


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
                #if int(row['cdr1.start']) > -1 and int(row['cdr1.end']) > -1 and int(row['cdr2.start']) > -1 and int(row['cdr2.end']) > -1 and int(row['cdr2.5.start']) > -1 and int(row['cdr2.5.end']) > -1 and ('TRAV' in row['id'] or 'TRBV' in row['id']):
                #    if row['sequence'][row['cdr1.start']:row['cdr1.end']] != test_cdrs[row['#species']][row['id']]['cdr1aa']:
                #        print('cdr1', row['#species'], row['sequence'][row['cdr1.start']:row['cdr1.end']], test_cdrs[row['#species']][row['id']]['cdr1aa'], row['id'])
                #    else:
                #        print('yay1')
                #    if row['sequence'][row['cdr2.start']:row['cdr2.end']] != test_cdrs[row['#species']][row['id']]['cdr2aa']:
                #       print('cdr2', row['#species'], row['sequence'][row['cdr2.start']:row['cdr2.end']], test_cdrs[row['#species']][row['id']]['cdr2aa'], row['id'])
                #    else:
                #        print('yay2')
                #    if row['sequence'][row['cdr2.5.start']:row['cdr2.5.end']] != test_cdrs[row['#species']][row['id']]['cdr2.5aa']:
                #        print('cdr2.5', row['#species'], row['sequence'][row['cdr2.5.start']:row['cdr2.5.end']], test_cdrs[row['#species']][row['id']]['cdr2.5aa'], row['id'])
                #    else:
                #        print('yay2.5')
                writer.writerow({key: row[key] for key in out_columns})
write_aa_sequences()


def get_tcr_genes(inpfile='tcr_aa_sequences.txt'):
    tcr_dict = {}
    with open(inpfile, 'r') as inp:
        reader = csv.DictReader(inp, delimiter='\t')
        for row in reader:
            if row['#species'] not in tcr_dict:
                tcr_dict[row['#species']] = {}
            tcr_dict[row['#species']][row['id']] = row['sequence']
    return(tcr_dict)


tcr_dict = get_tcr_genes()

def get_tcr(v_gene, j_gene, cdr3, species, tcr_dict=tcr_dict):
    if j_gene == '':
        return(tcr_dict[species][v_gene][:-1]+cdr3)
    else:
        return(tcr_dict[species][v_gene][:-1]+cdr3+tcr_dict[species][j_gene][1:])

