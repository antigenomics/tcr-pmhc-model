import urllib.request
import urllib.parse
import argparse
import os

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will parse aminoacid sequences for TR genes from various species from imgt")

parser.add_argument("-o", nargs=1, type=str, default="../imgt_work/imgt_all_aminoacid.fasta",#_nucleotide.fasta",
                    help="output file path.")
parser.add_argument("--seqtype", nargs=1, type=str, default="aminoacid", help="nucleotide or aminoacid sequences?")

args = parser.parse_args()
if type(args.o) is list:
    args.o = args.o[0]
output_file = args.o
if type(args.seqtype) is list:
    args.seqtype = args.seqtype[0]
seqtype = {"aminoacid": "7.3", "nucleotide": "7.1"}[args.seqtype]


outputinfo = []
tr_species = {'Homo+sapiens': 'Homo sapiens',
              'Mus+musculus': 'Mus musculus',
              'Bos+taurus': 'Bos taurus',
              'Canis+lupus+familiaris': 'Canis lupus familiaris',
              'Macaca+mulatta': 'Macaca mulatta',
              'Mus+minutoides': 'Mus minutoides',
              'Mus+pahari': 'Mus pahari',
              'Mus+spretus': 'Mus spretus',
              'Oncorhynchus+mykiss': 'Oncorhynchus mykiss',
              'Ovis+aries': 'Ovis aries'}

tr_group = ['TRAV', 'TRAJ', 'TRAC',
            'TRBV', 'TRBD', 'TRBJ', 'TRBC',
            'TRGV', 'TRGJ', 'TRGC',
            'TRDV', 'TRDD', 'TRDJ', 'TRDC']

#http://www.imgt.org/genedb/GENElect?query=7.3+TRAV&species=Homo+sapiens
#http://www.imgt.org/genedb/GENElect?query=7.3+Group&species=Species
#The FASTA header contains 15 fields separated by '|':
#
#1. IMGT/LIGM-DB accession number(s)
#2. IMGT gene and allele name
#3. species
#4. IMGT allele functionality
#5. exon(s), region name(s), or extracted label(s)
#6. start and end positions in the IMGT/LIGM-DB accession number(s)
#7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
#8. codon start, or 'NR' (not relevant) for non coding labels
#9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
#10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
#11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
#12. number of amino acids (AA): this field indicates that the sequence is in amino acids
#13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
#14. partial (if it is)
#15. reverse complementary (if it is)

spiti = 1
print(args.seqtype, seqtype)
for sp in tr_species:
    print(sp, str(spiti)+'/'+str(len(tr_species)))
    spiti += 1
    trgriti = 1
    for trgr in tr_group:
        print(trgr, str(trgriti)+'/'+str(len(tr_group)))
        trgriti += 1
        inpurl = urllib.request.urlopen('http://www.imgt.org/genedb/GENElect?query={}+{}&species={}'.format(seqtype, trgr, sp))
        inpurlheader = inpurl.headers.get_content_charset()
        for i in inpurl:
            if str(i).startswith('b\'<b>Number'): #part of the site, where it says how many entries are on the list
                for i in inpurl:
                    if str(i).startswith('b\'<pre>'): #next line will be fasta
                        for i in inpurl:
                            if str(i).startswith('b\'\\r\\n'): #end of fasta
                                break
                            else:
                                info = i.decode(inpurlheader)
                                if info.startswith('>'):
                                    info = info.split('|')
                                    info[2] = sp#tr_species[sp]
                                    info = '|'.join(info)
                                outputinfo.append(info)
                        break
                break
    print('\n')

with open(output_file, 'w') as out:
    out.writelines(outputinfo)