#!/usr/bin/env python3

import pytest
import importlib

from functools import partial

def _test(module_name=None, error_type=None, module_test_callback=None, assert_callback=None, package_name=None):
  if module_name is None or \
    error_type is None or \
    module_test_callback is None or \
    assert_callback is None or \
    package_name is None:
    raise ValueError('Insufficient arguments supplied to helper._test()')

  with pytest.raises(error_type) as e:
    full_module_name = '.'.join([package_name, module_name])
    module_test_callback(importlib.import_module(full_module_name))
  assert assert_callback(e)

def _test_error_str(module_name, error_type=None, error_str_search=None, package_name=None):
  _test(
    module_name=module_name, 
    error_type=error_type, 
    module_test_callback=lambda mod: mod.main(),
    assert_callback=lambda re: error_str_search in str(re.value), 
    package_name=package_name
  )

_test_script_main = partial(
    _test,
    error_type=SystemExit, 
    module_test_callback=lambda mod: mod.main(['--help']),
    assert_callback=lambda se: se.value.code == 0 # --help should return 0
  )

_test_non_executable = partial(
    _test_error_str,
    error_type=RuntimeError, 
    error_str_search='executed'
  )