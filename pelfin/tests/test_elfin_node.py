import pytest

def test_elfin_node():
  from pelfin.elfin_node import ElfinNode

  # Error and warning tests
  with pytest.raises(KeyError, match='name'):
    ElfinNode(id=-1)

  with pytest.raises(ValueError, match='Node ID.*'):
    ElfinNode(
      id=-1,
      name='Whatever'
    )

  with pytest.warns(UserWarning, match='.*trim vector.*'):
    ElfinNode(
      id=1, 
      name='Peter', 
      trim=[False, False],
      cap=[False, False],
      cterm_node_id=[]
    ) is not None

  with pytest.raises(ValueError, match='.*trimmed end.*'):
    ElfinNode(
      id=1, 
      name='Peter', 
      trim=[False, True],
      cap=[False, True],
      cterm_node_id=[]
    ) is not None