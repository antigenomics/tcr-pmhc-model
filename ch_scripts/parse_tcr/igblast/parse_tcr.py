import sys
import os
import argparse

import pandas as pd
from Bio import SeqIO

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(
    description="This script will find CDR positions in TCR sequences. "
                "Write IGDATA path to igblast , if you run it from different directory")

parser.add_argument(
    "--inpfile", nargs=1, type=str, default="../test_tcr_human.txt",#_nucleotide.fasta",
                    help="input file path.")
parser.add_argument(
    "--outfile", nargs=1, type=str, required=False, help="path to output data"
)
parser.add_argument(
    "--seqtype", nargs=1, type=str, default="aminoacid", help="nucleotide or aminoacid sequences?")
parser.add_argument(
    "--specie", nargs=1, type=str, default="HomoSapiens", help="Specie name in CamelCase")


igblast_db_path = os.path.dirname(os.path.realpath(__file__))
print("Running parse_tcr.py in {}\n".format(igblast_db_path))
args = parser.parse_args()

if type(args.inpfile) is list:
    args.inpfile = args.inpfile[0]
inpfile = args.inpfile

if type(args.seqtype) is list:
    args.seqtype = args.seqtype[0]
seq_type = args.seqtype

if type(args.specie) is list:
    args.specie = args.specie[0]
specie = args.specie

if args.outfile is not None:
    if type(args.outfile) is list:
        args.outfile = args.outfile[0]
    outfile = args.outfile
else:
    outfile = ".".join(inpfile.split(".")[:-1])+"_parsed.txt"
    

tmpfile_v = ".".join(inpfile.split(".")[:-1])+"_v_blasted.txt"
tmpfile_j = ".".join(inpfile.split(".")[:-1])+"_j_blasted.txt"

segments = ("Variable", "Diversity", "Joining")

workdir = os.path.join(igblast_db_path, "igblast_db", specie.lower(), "{}_{}_{}_sequences.fasta".format(seq_type, specie, "{}"))

v_db = workdir.format("Variable")
d_db = workdir.format("Diversity")
j_db = workdir.format("Joining")


if seq_type == "nucleotide":
    igblastcmd = "/usr/local/ncbi/igblast/bin/igblastn"
    blastcmd = "blastn"
    cmd = "-germline_db_J {} -germline_db_D {} ".format(j_db, d_db)
elif seq_type == "aminoacid":
    igblastcmd = "/usr/local/ncbi/igblast/bin/igblastp"
    blastcmd = "blastp"
    cmd = ""
else:
    sys.exit(1)


blastcmd = "{} -db {} -query {} -out {}".format(blastcmd, j_db, inpfile, tmpfile_j)
igblastcmd = "{} -germline_db_V {} {}-query {} -out {}".format(igblastcmd, v_db, cmd, inpfile, tmpfile_v)

print("Running BLAST:\n{}".format(blastcmd))
os.system(blastcmd)
print("Done\n")
print("Running IgBLAST:\n{}".format(igblastcmd))
os.system(igblastcmd)
print("Done\n")

segment_dict = {"FR1-IMGT": "FR1", "FR2-IMGT": "FR2", "CDR1-IMGT": "CDR1", "CDR2-IMGT":"CDR2"}
outinfo = []
with open(tmpfile_v, "r") as inp:
    for line in inp:
        if line.startswith("Query="):
            query = line[6:].strip()

            for line in inp:
                if line.startswith("\n"):
                    break
                query += line.strip()

            for line in inp:
                if line.startswith("Sequences producing significant alignments: "):
                    inp.readline()
                    v = inp.readline().strip().split()[0]
                if line.startswith(("FR1", "CDR1", "FR2", "CDR2", "FR3")):
                    info = line.strip().split()
                    start_pos = int(info[1])
                    end_pos = int(info[2])
                    if line.startswith("FR3"):
                        outinfo.append([query, "FR3", start_pos - 1, end_pos-1, v, ""])
                        outinfo.append([query, "CDR3", end_pos - 1, 0, v, ""])
                        break
                    elif info[0] in segment_dict:
                        outinfo.append([query, segment_dict[info[0]], start_pos - 1, end_pos, v, ""])


outdf = pd.DataFrame(outinfo, columns=["query", "region", "start", "end", "v.gene", "j.gene"])

with open(tmpfile_j, "r") as inp:
    for line in inp:
        if line.startswith("Query="):
            query = line[6:].strip()
            for line in inp:
                if line.startswith("\n"):
                    break
                query += line.strip()
            for line in inp:
                if line.startswith("***** No hits found *****"):
                    outdf.loc[(outdf["query"] == query) & (outdf["region"] == "CDR3"), "end"] = -1
                    outdf.loc[(outdf["query"] == query), "j.gene"] = -1
                    break
                if line.startswith(">"):
                    j = line.strip().split()[0].replace(">", "")
                    outdf.loc[(outdf["query"] == query), "j.gene"] = j
                if line.startswith("Query "):
                    info = line.split()
                    outdf.loc[(outdf["query"] == query) & (outdf["region"] == "CDR3"), "end"] = int(info[1])
                    break

#outdf["end"] - outdf["start"]


def slice_sequence(df, fullseqcol, seqcol, startcol, endcol):
    def slice_seq(fullseq, start, end):
        return fullseq[start:end]

    df.loc[:, seqcol] = df.apply(lambda col: slice_seq(col[fullseqcol], col[startcol], col[endcol]), axis=1)
    return df[seqcol]


with open(inpfile, "r") as inp:
    for record in SeqIO.parse(inp, "fasta"):
        outdf.loc[(outdf["query"] == record.description), "fullseq"] = str(record.seq)


outdf["sequence"] = slice_sequence(outdf, "fullseq", "sequence", "start", "end")
outdf["length"] = outdf["sequence"].str.len()

outdf.to_csv(outfile, columns=["query", "region", "start", "end", "length", "sequence", "v.gene", "j.gene"], sep="\t", index=None)


#os.remove(tmpfile_j)
#os.remove(tmpfile_v)
