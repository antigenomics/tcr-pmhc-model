import os

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


c_domains = pd.read_table("c_domain_PDB_info.txt")
chains = {"TRBC": "TCRB", "TRAC": "TCRA"}

specdict = {"MusMusculus": "MOUSE", "HomoSapiens": "HUMAN"}
for spec in specdict:
    sequences = []
    with open("sequence_db/{}/aminoacid_{}_Constant_sequences.fasta".format(spec.lower(), spec), 'r') as inp:
        for record in SeqIO.parse(inp, format="fasta"):
            chain = record.id[:4]
            longest_sequence = c_domains[(c_domains["specie"] == spec) &
                                         (c_domains["chain"] == chains[chain])].reset_index(drop=True)["sequence"][0]
            alignment = pairwise2.align.globaldx(longest_sequence, record.seq, matlist.blosum62)[0]
            seq1 = alignment[0]
            seq2 = alignment[1]
            for i in range(1, len(seq1)+1):
                if seq1[-i] != '-':
                    seq2 = seq2[:-i+1].replace("-", "")
                    break
            sequences.append(SeqRecord(Seq(seq2, IUPAC.protein), id=record.id, description=""))
    with open("sequence_db/{}/aminoacid_{}_Constant_sequences_splitted_for_pdb.fasta".format(spec.lower(), spec),
              'w') as out:
        SeqIO.write(sequences, out, format="fasta")




