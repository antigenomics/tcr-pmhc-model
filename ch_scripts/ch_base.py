"""
Some lines of code were imported from tcr-dist package with little modifications (https://github.com/phbradley/tcr-dist)

MIT License

Copyright (c) 2017 Philip Harlan Bradley and Jeremy Chase Crawford

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os
import sys

class Found(Exception): pass

orgtoorg = {'HomoSapiens':'human', 'MusMusculus':'mouse'}
org_to_vdjdb_org = {'human':'HomoSapiens', 'Human':'HomoSapiens', 'homosapiens':'HomoSapiens', 'HomoSapiens':'HomoSapiens',
                    'mouse':'MusMusculus', 'Mouse':'MusMusculus', 'musmusculus':'MusMusculus', 'MusMusculus':'MusMusculus'}


#http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
anchors = {'CDR1': [26, 38], 'CDR2': [55, 65], 'CDR2.5': [80, 86]} #'CDR3': [104, 116]
mouse_anchors = {'CDR1': [27, 39], 'CDR2': [56, 66], 'CDR2.5': [81, 88]}


def crdir(path):
    if os.path.exists(path):
        return 0
    else:
        os.mkdir(path)

def get_filename_no_ext(input_file):
    return input_file.rsplit('.', 1)[0]

def parse_tsv_line(tsvline, headinfo):
    if type(tsvline) != list:
        wtsvline = tsvline.split('\t')
    else:
        wtsvline = tsvline
    wtsvline[-1] = wtsvline[-1].strip()
    assert(len(wtsvline)==len(headinfo))
    outinfo = {}
    for tag,val in zip(headinfo, wtsvline):
        outinfo[tag] = val
    return outinfo

def parse_listtsv_line(tsvline, headinfo):
    outinfo = {}
    for tag,val in zip(headinfo, tsvline):
        outinfo[tag] = val
    return outinfo

def make_tsv_line(ininfo, headinfo, empty_string_replacement=''):
    l = []
    for tag in headinfo:
        val = ininfo[tag]
        if type(val) is str:
            if empty_string_replacement and not val:
                l.append(empty_string_replacement)
            else:
                l.append(val)
        else:
            l.append(str(val))
    return '\t'.join(l)


def ch_check(smth):
    print(smth)
    sys.exit()


BASE_SCRIPTDIR = os.path.dirname(os.path.abspath('ch_base.py'))
BASE_HOMEDIR = BASE_SCRIPTDIR[:BASE_SCRIPTDIR.rfind('AdaptiveImm')+len('AdaptiveImm')+1]