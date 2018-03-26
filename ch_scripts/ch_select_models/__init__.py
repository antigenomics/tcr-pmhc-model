import os, sys

init_dirname = os.path.dirname(__file__)

if __name__ == '__main__' and __package__ is None:
    print('Working in {}'.format(os.path.basename(init_dirname)))

def get_package():
    print('Initializing package: {}'.format(os.path.basename(init_dirname)))
    if os.path.split(os.path.split(init_dirname)[0])[0] not in sys.path:
        sys.path.insert(1, os.path.split(os.path.split(init_dirname)[0])[0])
    if os.path.split(init_dirname)[0] not in sys.path:
        sys.path.insert(1, os.path.split(init_dirname)[0])
    return 'ch_scripts.{}'.format((os.path.basename(init_dirname)))

if os.path.dirname(__file__) not in sys.path:
    sys.path.insert(1, init_dirname)