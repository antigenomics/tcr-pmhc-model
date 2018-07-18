from Bio import SeqIO
import pandas as pd

specdict = {"Homo+sapiens": "HomoSapiens",
            "Mus+musculus": "MusMusculus",
            "Macaca+mulatta": "MacacaMulatta"}


# for Constant regions:
# imgt.org/IMGTrepertoire/Proteins/alleles/mouse/TRA/TRAC/Mu_TRAC.html


with open("imgt_all_nucleotide.fasta", "r") as inp:
    dictinfo = {}
    for record in SeqIO.parse(inp, "fasta"):
        info = record.description.split("|")
        gene = info[1]
        spec = info[2]
        fnucleotide = info[8]
        if (gene.startswith("TRAC") or gene.startswith("TRBC") or gene.startswith("TRDC") or gene.startswith("TRGC")) \
                and spec in specdict and "untranslated" not in info[4]:
            if fnucleotide != "+1":
                print(gene, spec, info)
                print(fnucleotide)
            name = "{}\t{}".format(specdict[spec], gene)
            if name not in dictinfo:
                dictinfo[name] = ""
            dictinfo[name] += str(record.seq).replace(".", "").upper()
    outinfo = []
    for name in dictinfo:
        outinfo.append([*name.split("\t"), dictinfo[name]])
    pd.DataFrame(outinfo, columns=["species", "c.id", "c.seq"]).to_csv("constant.nucleotide.txt", sep="\t", index=None)

