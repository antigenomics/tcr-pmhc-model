import sys
import os
import argparse


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
    outfile = ".".join(inpfile.split(".")[:-1])+"_igblast_out.txt"

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


blastcmd = "{} -db {} -query {} -out {}".format(blastcmd, j_db, inpfile, "")

cmd = "{} -germline_db_V {} {}-query {} -out {}".format(igblastcmd, v_db, cmd, inpfile, outfile)
print(cmd)
os.system(cmd)

