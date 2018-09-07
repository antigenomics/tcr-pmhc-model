import os
import sys
import argparse


#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(
    description="This script will create blast db from sequences.")

parser.add_argument(
    "--inpdb", nargs=1, type=str, default="sequence_db",#_nucleotide.fasta",
                    help="input file path.")
parser.add_argument(
    "--outdb", nargs=1, type=str, default="igblast/igblast_db", help="path to output data"
)
parser.add_argument(
    "--seqtype", nargs=1, type=str, default="aminoacid", help="nucleotide or aminoacid sequences?")

args = parser.parse_args()
if type(args.inpdb) is list:
    args.inpdb = args.inpdb[0]
sequencedb = args.inpdb
if type(args.outdb) is list:
    args.outdb = args.outdb[0]
outdb = args.outdb
if type(args.seqtype) is list:
    args.seqtype = args.seqtype[0]
seq_type = args.seqtype


def crdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

crdir(outdb)

species = ("HomoSapiens", "MusMusculus")

for seq_type in ("aminoacid", "nucleotide"):
    for spec in species:
        crdir(os.path.join(outdb, spec.lower()))
        for segm in ("Variable", "Diversity", "Joining", "Constant"):
            workfile = os.path.join(spec.lower(), "{}_{}_{}_sequences.fasta".format(seq_type, spec, segm))
            inpfile = os.path.join(sequencedb, workfile)
            outfile = os.path.join(outdb, workfile)

            dbtype = {"nucleotide":"nucl", "aminoacid":"prot"}[seq_type]

            cmd = "makeblastdb -parse_seqids -dbtype {} -in {} -out {}".format(dbtype, inpfile, outfile)
            os.system(cmd)