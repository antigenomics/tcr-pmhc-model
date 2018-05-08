import os
import sys

class Found(Exception): pass

def crdir(path):
    if os.path.exists(path):
        return 0
    else:
        os.mkdir(path)


def ch_check(smth):
    print(smth)
    sys.exit()


BASE_SCRIPTDIR = os.path.dirname(os.path.abspath('ch_base.py'))
BASE_HOMEDIR = BASE_SCRIPTDIR[:BASE_SCRIPTDIR.rfind('AdaptiveImm')+len('AdaptiveImm')+1]