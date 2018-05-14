#!/usr/bin/env python

import ElfinUtils
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 


def main():
	# Note: inf length ~22
	machines = [
		'IvyBridge\nZoo\n24-core\n2.7Ghz',
		'KNL\nZoo\n64-core\n1.3Ghz',
		'JetsonTX1\nZoo\n4-core\n1.9Ghz',
		'Core i7\nMacbookPro\n4-core\n2.6Ghz',
		'SandyBridge\nBC3\n16-core\n2.6Ghz',
		'SandyBridge\nSwan\n8-core\n2.6Ghz'
	]

	# Intel 2s12c1t IvyBridge Xeon E5-2697 v2 2.7Ghz
	ivy1t = 26240.0
	ivy24t = 1344.0
	ivyCores = 24
	ivyInfTime = 38633 #score 8255

	# Intel 1s64c4t Xeon Phi Kights Landing 1.3Ghz
	knl1t = 128642.0
	knl64t = 2479.0
	knlCores = 64
	knlInfTime = 58713 #score 8255

	# nVidia 4s1c1t Jetson TX1 1.9Ghz
	jet1t = 59776.0
	jet4t = 14929.0
	jetCores = 4
	jetInfTime = 6*60000+56401 #score 8682

	# MBP Intel 1s4c2t IvyBridge Core i7 2.6Ghz
	mbp1t = 25204.0
	mbp4t = 8000.0
	mbpCores = 4
	mbpInfTime = 4*60000+58061 #score 8696

	# BC3 Intel 2s8c1t SandyBridge 2.6Ghz
	bc1t = 27211.0
	bc16t = 1905.0
	bcCores = 16
	bcInfTime = 48027 #score 9215

	# Swan Intel 1s8c2t SandyBridge 2.6Ghz
	swa1t = 28566.0
	swa8t = 4313.0
	swaCores = 8
	swaInfTime = 1*60000+47078 #score 10433


	speedUps = [
		ivy1t/ivy24t,
		knl1t/knl64t,
		jet1t/jet4t,
		mbp1t/mbp4t,
		bc1t/bc16t,
		swa1t/swa8t
	]

	cores = [
		ivyCores,
		knlCores,
		jetCores,
		mbpCores,
		bcCores,
		swaCores
	]

	times = [ t / 1000.0 for t in [
		ivyInfTime,
		knlInfTime,
		jetInfTime,
		mbpInfTime,
		bcInfTime,
		swaInfTime
	]]

	speedUpsPerCore = [s/c for (s,c) in zip(speedUps, cores)]

	N = len(speedUpsPerCore)
	ind = np.arange(N)
	width = 0.7
	hOffsetPerc = 0.01
	fontSize = 16
	axesArea = [0.1, 0.15, .8, .8]
	figSize = (24, 9)
	transparent = True

	fig1 = plt.figure(figsize=figSize)
	ax1 = plt.axes(axesArea, frameon=True)

	# Beautiful color switching code thanks to user1839053's SO post at
	# http://stackoverflow.com/questions/4971269/how-to-pick-a-new-color-for-each-plotted-line-within-a-figure-in-matplotlib
	bars = ax1.bar(ind, speedUpsPerCore, width)
	color=iter(cm.rainbow(np.linspace(0,1,N)))
	for b in bars:
		b.set_color(next(color))

	ax1.set_ylabel('Normalised OpenMP Speedup', fontsize=fontSize+5)

	plt.xticks(ind, machines)
	ax1.tick_params(axis='both', which='minor', labelsize=fontSize)
	ax1.tick_params(axis='both', which='major', labelsize=fontSize)

	rects1 = ax1.patches
	labels1 = [('(' + str(round(s, 1)) + 'x)') for s in speedUps]
	hOffset = max(speedUpsPerCore) * hOffsetPerc;
	for rect, label in zip(rects1, labels1):
	    height = rect.get_height() + hOffset
	    ax1.text(rect.get_x() + rect.get_width()/2, height, label, ha='center', va='bottom', fontsize=fontSize)

	fig1.show()
	plt.savefig('OpenMP.png', transparent=transparent)

	# Hijack

	machines = [
		'Intel i7\n3720\nMacbookPro\n4-core\n2.6Ghz',
		'Intel i7\n4790\nZoostorm\n4-core\n2.6Ghz',
		'Intel Xeon\nE5-2697 v2\n(Zoo)\n24-core\n2.7Ghz',
		'Intel Xeon\nPhi 7210\n(Zoo)\n64-core\n1.3Ghz',
		'Intel Xeon\nE5-2670\n(BC3)\n16-core\n2.6Ghz',
		'Intel Xeon\nE5-2670\n(Swan)\n8-core\n2.6Ghz',
		'Intel Xeon\nE5-2640\n(Bluegem)\n16-core\n2.6Ghz',
		'Nvidia Tesla\nK40m\n(Zoo)\n4.3TFlops',
		'Nvidia GTX\n980 Ti\n(Zoo)\n5.6TFlops',
		'Nvidia GTX\n1080 Ti\n(Zoo)\n10.6TFlops'
	]

	times = [
		1110.3,
		779.2,
		209.9,
		269.7,
		261.5,
		448.2,
		293.8,
		234.8,
		227.1,
		167.7
	]

	N = len(times)
	ind = np.arange(N)

	# Plot time
	fig2 = plt.figure(figsize=figSize)
	ax2 = plt.axes(axesArea, frameon=True)

	bars = ax2.bar(ind, times, width)
	color=iter(cm.rainbow(np.linspace(0,1,N)))
	for b in bars:
		b.set_color(next(color))

	ax2.set_ylabel('Time-to-solution (s)', fontsize=fontSize+5)

	plt.xticks(ind, machines)
	ax2.tick_params(axis='both', which='minor', labelsize=fontSize)
	ax2.tick_params(axis='both', which='major', labelsize=fontSize)

	rects2 = ax2.patches
	labels2 = [round(t, 1) for t in times]
	hOffset = max(times) * hOffsetPerc;
	for rect, label in zip(rects2, labels2):
	    height = rect.get_height() + hOffset
	    ax2.text(rect.get_x() + rect.get_width()/2, height, label, ha='center', va='bottom', fontsize=fontSize)

	fig2.show()
	plt.savefig('TimeToSol.png', transparent=transparent)

	input()

if __name__ == '__main__':
	ElfinUtils.safeExec(main)