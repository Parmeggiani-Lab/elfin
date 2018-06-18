#!/usr/bin/env python

import argparse, sys
# from ElfinUtils import *

def main():
	ap = argparse.ArgumentParser(description='Template Elfin Python script')
	ap.add_argument('input') # Absence of dash denotes mandatory argument
	args = ap.parse_args()

if __name__ == '__main__':
	main()