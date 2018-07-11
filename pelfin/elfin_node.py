import numpy as np
import warnings

class ElfinNode():
  '''
  A single module instance and stores info about connectivity
  '''
  def __init__(
    self, 
    **kwargs
  ):
    self.id = kwargs.pop('id')
    self.name = kwargs.pop('name')
    self.trim = kwargs.pop('trim', (False, False))
    self.cap = kwargs.pop('cap', None)
    self.cterm_node_id = kwargs.pop('cterm_node_id', -1)
    self.rot = kwargs.pop('rot', [[1,0,0],[0,1,0],[0,0,1]])
    self.tran = kwargs.pop('tran', [0,0,0])

    # Default is to cap the end that is not trimmed
    if self.cap == None:
      self.cap = (not self.trim[0], not self.trim[1])

    # Error checking
    if self.id < 0:
      raise ValueError('Node ID should never be negative: id={}'.format(self.id))

    if len(self.trim) != 2:
      raise ValueError('ElfinNode trim vector length != 2: trim={}'.format(self.trim))

    if not self.trim[0] and not self.trim[1]:
      warnings.warn('ElfinNode trim vector both ends are NOT trimmed. '
      'This should only happen if the chain has one single node, which '
      'is not thought to be common. Proceed only if this is deliberate.')

    for i in [0, 1]:
      if self.trim[i] and self.cap[i]:
        raise ValueError('Cannot cap a trimmed end[{}]: name={}, id={}'
          .format(i, self.name, self.id))

  def transform(self, rot, tran):
      self.rot = (np.dot(self.rot, rot)).tolist()
      self.tran = (np.dot(self.tran, rot) + tran).tolist()

def main():
  raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
  main()