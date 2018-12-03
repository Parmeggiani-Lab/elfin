class ElfinGraph(object):
    """
    A network of nodes that are connected either by doubles or through hubs.
    Might be multi-chain.
    """
    def __init__(self, name='', nodes=[]):
        self.name = name
        self.nodes = nodes

    def __repr__(self):
        return 'ElfinGraph: {{\n{}\n}}\n'.format(
            '\n'.join((repr(n) for n in self.nodes)))    
        
    def transform(self, rot, tran):
        for n in self.nodes:
            n.transform(rot, tran)


def main():
    """main"""
    raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
    main()