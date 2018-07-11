import pytest
import importlib

from functools import partial

from tests import helper

test_package_name = 'pelfin'
script_main_test = partial(helper._test_script_main, package_name=test_package_name)
non_executable_test = partial(helper._test_non_executable, package_name=test_package_name)

def test_mains():
  # Test python code in order of workflow
  script_main_test('template')

  non_executable_test('utilities')
  non_executable_test('pdb_utilities')
  non_executable_test('elfin_graph')
  non_executable_test('elfin_node')
  non_executable_test('kabsch')

  script_main_test('preprocess')
  script_main_test('hubinfo_convert')
  script_main_test('v1_design_convert')
  script_main_test('dbgen')
  script_main_test('stitch')

  script_main_test('job_dispatcher')
  script_main_test('rmsd')
  script_main_test('stat_xdb')