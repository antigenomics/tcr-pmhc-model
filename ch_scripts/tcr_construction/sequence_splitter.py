import itertools
import argparse
import os
import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(
    description="This script will split parsed aminoacid "
                "sequences by types of their regions.")

parser.add_argument(
    "-i", nargs=1, type=str,#_nucleotide.fasta",
                    help="input file path.", required=False)
parser.add_argument(
    "-outdir", nargs=1, type=str, default="sequence_db", help="path to output data"
)
parser.add_argument(
    "--seqtype", nargs=1, type=str, default="aminoacid", help="nucleotide or aminoacid sequences?")


args = parser.parse_args()

if type(args.outdir) is list:
    args.outdir = args.outdir[0]
outdir = args.outdir
if type(args.seqtype) is list:
    args.seqtype = args.seqtype[0]
seq_type = args.seqtype


def crdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

if args.i is not None:
    if type(args.i) is list:
        args.i = args.i[0]
    sequences = pd.read_table(args.i)

else:
    if seq_type == "nucleotide":
        inpfile_vj = "tcr_genes.txt"
        inpfile_c = "constant.nucleotide.txt"
    elif seq_type == "aminoacid":
        inpfile_vj = "tcr_aa_sequences.txt"
        inpfile_c = "constant.txt"
    else:
        sys.exit(1)

    sequences_vj = pd.read_table(inpfile_vj)
    sequences_c = pd.read_table(inpfile_c)
    sequences_c.columns = ["#species", "id", "sequence"]
    sequences_c["segment"] = "Constant"
    sequences = sequences_vj.append(sequences_c)

if seq_type == "nucleotide":
    seq_record_type = "generic_dna"
elif seq_type == "aminoacid":
    seq_record_type = "generic_protein"
else:
    seq_record_type = seq_type


def SeqRecord_maker(df, seq_record_col, seq_col, seq_type, id_col, description=""):
    if seq_type == "nucleotide":
        seq_record_type = "generic_dna"
    elif seq_type == "aminoacid":
        seq_record_type = "generic_protein"
    else:
        seq_record_type = seq_type


    def get_SeqRecord(sequence, id):
        return SeqRecord(Seq(sequence, seq_record_type), id=id, description=description)

    df[seq_record_col] = df.apply(lambda col: get_SeqRecord(col[seq_col], col[id_col]), axis=1)
    return df


SeqRecord_maker(sequences, "SeqRec", "sequence", seq_record_type, "id")

crdir(outdir)
for specie, specseq in sequences.groupby("#species"):
    specdir = os.path.join(outdir, specie.lower())
    crdir(specdir)
    for segment, segmseq in specseq.groupby("segment"):
        print(list(segmseq["SeqRec"]))
        with open(os.path.join(specdir, "{}_{}_{}_sequences.fasta".format(seq_type, specie, segment)), "w") as out:
            SeqIO.write(list(segmseq["SeqRec"]), out, "fasta")