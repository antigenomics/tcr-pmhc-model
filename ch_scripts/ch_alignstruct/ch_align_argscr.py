#initiate package

if __name__ == '__main__' and __package__ is None:
    import sys
    sys.path.insert(0, '../..')
    sys.path.insert(0, '..')
    print(sys.path)
    print('"__package__" is None.')
    __package__ = 'ch_scripts.ch_alignstruct'


from . import ch_get_superimpose_structures