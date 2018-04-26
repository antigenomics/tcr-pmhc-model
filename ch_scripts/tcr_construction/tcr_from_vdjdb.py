import os
import csv
from gene_to_protein import tcr_dict, get_tcr

BASE_SCRIPTDIR = os.path.dirname(os.path.abspath('tcr_from_vdjdb.py'))
BASE_HOMEDIR = BASE_SCRIPTDIR[:BASE_SCRIPTDIR.rfind('AdaptiveImm')+len('AdaptiveImm')+1]
filter_spec = 'HomoSapiens'
tcrs = set()

out_columns = ['complex.id', 'gene', 'cdr3', 'v.segm', 'j.segm', 
           'species', 'mhc.a', 'mhc.b', 'mhc.class', 
           'antigen.epitope' ,'tcr.sequence']
with open('tcrs_{}_from_vdjdb.txt'.format(filter_spec), 'w') as out:
    writer = csv.DictWriter(out, fieldnames=out_columns, delimiter='\t')
    writer.writeheader()
    with open(os.path.join(BASE_HOMEDIR, 'vdjdb-2018-02-03/vdjdb.txt'), 'r') as inp:
        reader = csv.DictReader(inp, delimiter='\t')
        for row in reader:
            if row['species'] == filter_spec and row['v.segm'].endswith('1') and row['j.segm'].endswith('1') and row['cdr3'] != '':
                print([row['v.segm'], row['j.segm'], row['cdr3']])
                row['tcr.sequence'] = get_tcr(row['v.segm'], row['j.segm'], row['cdr3'], row['species'], tcr_dict)
                writer.writerow({key: row[key] for key in out_columns})
