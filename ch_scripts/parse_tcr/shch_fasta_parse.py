import os
specdict = {"MusMusculus":"MOUSE", "HomoSapiens":"HUMAN"}
os.chdir("igblast")

for spec in specdict:
    for chain in ("TCRA", "TCRB"):
        inpfile = "../{}_PDB_{}.fasta".format(specdict[spec], chain)
        cmd = "source activate py36 ;python3 parse_tcr.py --inpfile {} --seqtype {} --specie {}".format(inpfile,
                                                                                          "aminoacid",
                                                                                          spec)
        os.system(cmd)
        inpfile2 = "../CHIMERA_PDB_{}.fasta".format(chain)
        outfile2 = "../CHIMERA_{}_PDB_{}_parsed.txt".format(specdict[spec], chain)
        cmd2 = "source activate py36; python3 parse_tcr.py --inpfile {} --outfile {} --seqtype {} --specie {}".format(inpfile,
                                                                                                        outfile2,
                                                                                                        "aminoacid",
                                                                                                        spec)
        os.system(cmd2)

