import pytest
import importlib

def _test_script_main(script_name):
  with pytest.raises(SystemExit) as se:
    module_name = '.'.join(['pelfin', script_name])
    script = importlib.import_module(module_name)
    script.main(test_args=['--help'])
  assert se.value.code == 0

def _test_non_executable(non_exe_name):
  with pytest.raises(RuntimeError) as re:
    module_name = '.'.join(['pelfin', non_exe_name])
    non_exe = importlib.import_module(module_name)
    non_exe.main()
  assert 'executed' in str(re.value)

def test_mains():
  # Test python code in order of workflow
  _test_script_main('template')

  _test_non_executable('utilities')
  _test_non_executable('pdb_utilities')
  _test_non_executable('elfin_graph')
  _test_non_executable('elfin_node')
  _test_non_executable('kabsch')

  _test_script_main('preprocess')
  _test_script_main('hubinfo_convert')
  _test_script_main('v1_design_convert')
  _test_script_main('dbgen')
  _test_script_main('stitch')

  _test_script_main('job_dispatcher')
  _test_script_main('rmsd')
  _test_script_main('stat_xdb')