import sys
import os
import argparse
import pandas as pd

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(
    description="This search TCR sequence in .")

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
    outfile = ".".join(inpfile.split(".")[:-1])+"_parsed_out.txt"
    

tmpfile_v = ".".join(inpfile.split(".")[:-1])+"_v_blasted.txt"
tmpfile_j = ".".join(inpfile.split(".")[:-1])+"_j_blasted.txt"

segments = ("Variable", "Diversity", "Joining")

workdir = os.path.join("igblast_db", specie.lower(), "{}_{}_{}_sequences.fasta".format(seq_type, specie, "{}"))

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
os.system(blastcmd)
os.system(igblastcmd)


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
                if line.startswith("CDR1"):
                    info = line.strip().split()
                    outinfo.append([query, "CDR1", int(info[1])-1, int(info[2])])
                if line.startswith("CDR2"):
                    info = line.strip().split()
                    outinfo.append([query, "CDR2", int(info[1]) - 1, int(info[2])])
                if line.startswith("FR3"):
                    info = line.strip().split()
                    outinfo.append([query, "CDR3", int(info[2]) - 1, 0])
                    break


outdf = pd.DataFrame(outinfo, columns=["query", "CDRtype", "start", "end"])


with open(tmpfile_j, "r") as inp:
    for line in inp:
        if line.startswith("Query="):
            query = line[6:].strip()
            for line in inp:
                if line.startswith("\n"):
                    break
                query += line.strip()
            for line in inp:
                if line.startswith("Query "):
                    info = line.split()
                    outdf.loc[(outdf["query"] == query) & (outdf["CDRtype"] == "CDR3"), ("end")] = int(info[1])
                    break
outdf.to_csv(outfile, sep="\t", index=None)


os.remove(tmpfile_j)
os.remove(tmpfile_v)
