#!/usr/bin/env python

import ElfinUtils
import json
import argparse
import numpy as np
from decimal import Decimal

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def main():
	ap = argparse.ArgumentParser(description='Compute the number of combinations for a given MMC protein length');
	ap.add_argument('--xdbFile', default='./resources/xDB.json')
	ap.add_argument('--length', type=int, default=21)

	globals().update(vars(ap.parse_args()))

	print xdbFile
	with open(xdbFile, 'r') as file:
		xdb = json.load(file)
		dd = xdb['doublesData']

		singleNames = xdb['singlesData'].keys()
		dim = len(singleNames)

	dim = len(singleNames)
	adjMat = np.zeros([dim, dim])

	for pdk in dd.keys():
		i1 = singleNames.index(pdk)
		for pdkk in dd[pdk].keys():
			i2 = singleNames.index(pdkk)
			adjMat[i1][i2] = 1

	mmcYs = []
	for l in xrange(1, length):
		nCombs = Decimal(np.sum(np.linalg.matrix_power(adjMat, l)))

		mmcYs.append(Decimal(nCombs))
		# print 'L={}, NC={}'.format(l+1, nCombs)

 	# A typical repeat module is ~100 AA
	Xs = np.asarray(range(2, length + 1))
	Xs = [x*100 for x in Xs]
	mmcYs = np.asarray(mmcYs)
	aaYs = np.power(20.0, np.asarray(Xs, dtype=np.float64))

	fig, ax1 = plt.subplots()

	ax1.set_xlabel('Design Length/AA')
	# ax2 = ax1.twinx()

	ElfinUtils.pauseCode()
	ax1.plot(Xs[:len(aaYs)], aaYs, label='AA')
	ax1.set_ylabel('No. Combs (log scale)')
	ax1.set_yscale('log')
 	
	ax1.plot(Xs, mmcYs, label='MMC')

	# xTickIds = np.arange(3, len(Xs) + 1, 5)
	# xTickIds = np.insert(xTickIds, 0, 0)
	# plt.xticks(xTickIds+2, [Xs[xtId] for xtId in xTickIds])

	plt.legend()
	plt.show()

if __name__ == '__main__':
	main()