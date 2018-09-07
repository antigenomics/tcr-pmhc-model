import sys
import os
import argparse

import pandas as pd
from Bio import SeqIO

specdict = {"MusMusculus": "MOUSE", "HomoSapiens": "HUMAN"}

outinfo = []


for spec in specdict:
    for chain in ("TCRA", "TCRB"):
        c_db = os.path.join("igblast/igblast_db", spec.lower(), "{}_{}_{}_sequences.fasta".format("aminoacid", spec, "Constant"))
        seqfile = "{}_PDB_{}.fasta".format(specdict[spec], chain)
        outfile = "{}_PDB_{}_constant_tmp.txt".format(specdict[spec], chain)
        blastcmd = "blastp -db {} -query {} -out {}".format(c_db, seqfile, outfile)
        os.system(blastcmd)
        longest_c = [0, ""]
        with open(outfile, 'r') as inp:
            for line in inp:
                if line.startswith("Query="):
                    for line in inp:
                        if "Positives" in line:
                            info = line.strip().split()
                            alignlength = int(info[6].split("/")[1])
                            if alignlength > longest_c[0]:
                                seq = ""
                                for line in inp:
                                    if line.startswith("Sbjct "):
                                        seq += line.strip().split()[2].replace("-", "")
                                    elif line.startswith(("Lambda", ">")):
                                        longest_c = [alignlength, seq]
                                        break
                            break
        outinfo.append([spec, chain, longest_c[0], longest_c[1]])

outdf = pd.DataFrame(outinfo, columns=["specie", "chain", "length", "sequence"])
outdf.to_csv("c_domain_PDB_info.txt", sep="\t", index=None)
