import pytest
import importlib

from functools import partial

from tests import helper

test_package_name = 'pymol_extensions'
script_main_test = partial(helper._test_script_main, package_name=test_package_name)
non_executable_test = partial(helper._test_non_executable, package_name=test_package_name)

def test_mains():
  non_executable_test('extension_template')
  non_executable_test('transform_helper')
  non_executable_test('load_all_extensions')
  non_executable_test('draw_lines')
  non_executable_test('batch_convert')