import sys, os

curr_dir = os.path.dirname(__file__)
sys.path.append(os.path.join(curr_dir, os.pardir)) # for _test_utilities
sys.path.append(os.path.join(curr_dir, os.pardir, os.pardir, 'pelfin')) # for test scripts