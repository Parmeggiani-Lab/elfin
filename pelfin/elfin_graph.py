class ElfinGraph():
  """
  A network of nodes that are connected either by doubles or through hubs.
  Might be multi-chain.
  """
  def __init__(self, name='', nodes=[]):
    self.name = name
    self.nodes = nodes
    
  def transform(self, rot, tran):
    for n in self.nodes:
      n.transform(rot, tran)

def main():
  raise RuntimeError('This module should not be executed like a script')

if __name__ =='__main__': 
  main()