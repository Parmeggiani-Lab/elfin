import numpy as np
import warnings

class ElfinNode(object):
    """
    A single module instance and stores info about connectivity
    """
    def __init__(
        self, 
        **kwargs
    ):
        self.id = kwargs.pop('id')
        self.name = kwargs.pop('name')
        self.trim = kwargs.pop('trim', {'n': False, 'c': False})

        # Default is to cap the end that is not trimmed
        self.cap = kwargs.pop('cap', (not self.trim['n'], not self.trim['c']))
        self.cterm_node_id = kwargs.pop('cterm_node_id', -1)
        self.rot = kwargs.pop('rot', [[1,0,0],[0,1,0],[0,0,1]])
        self.tran = kwargs.pop('tran', [0,0,0])

        # Error checking
        if self.id < 0:
            raise ValueError('Bad ElfinNode id: {}'.format(self.id))

        if len(self.trim) != 2 or \
            'n' not in self.trim or 'c' not in self.trim:
            raise ValueError('Bad ElfinNode trimming flags: {}'.format(self.trim))

        if not self.trim['n'] and not self.trim['c']:
            warnings.warn('ElfinNode (id={}, name={}) '
                'trimming flags are FALSE for both n and c terms. '
                'This only happens if there is a single module network or chain. '
                'Proceed only if this is deliberate.'.format(self.id, self.name))

        for term in {'n', 'c'}:
            if self.trim[term] and self.cap[term]:
                raise ValueError('Cannot cap a trimmed end({}): name={}, id={}'
                    .format(term, self.name, self.id))

    def __repr__(self):
        return 'ElfinNode: ID={}, Name={}'.format(self.id, self.name)

    def transform(self, rot, tran):
            self.rot = (np.dot(self.rot, np.transpose(rot))).tolist()
            self.tran = (np.dot(self.tran, np.transpose(rot)) + tran).tolist()

def main():
    """main"""
    raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
    main()