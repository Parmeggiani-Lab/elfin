#!/usr/bin/env python

import ElfinUtils
import json
import argparse
from collections import OrderedDict
import GridSearchParams

def main():
	ap = argparse.ArgumentParser(description='Generate Grid Search configurations');
	ap.add_argument('--bmDir', default='./bm/l20/')
	ap.add_argument('--gsVersion', type=int, default=2)

	args = ap.parse_args()
	globals().update(vars(args))

	gsConfDir = './gsConfigsV{}/'.format(gsVersion)
	gsOutDir = './gsOutV{}/'.format(gsVersion)

	gsParams = ElfinUtils.Bunch(GridSearchParams.getGSParams(gsVersion))
	print 'Total runs needed: {}'.format(gsParams.nRuns)

	ElfinUtils.mkdir(gsConfDir)

	# Write all combinations of GA parameters to output
	configId = 0
	for cld in gsParams.chromoLenDevs:
		for gps in gsParams.gaPopSizes:
			for gi in gsParams.gaIters:
				for gsr in gsParams.gaSurviveRates:
					for gcr in gsParams.gaCrossRates:
						for (gpmr, glmr) in zip(gsParams.gaPointMutateRates, gsParams.gaLimbMutateRates):
							for bmName in gsParams.bmNames:
								outputName = 'gs_{}_{}'.format(configId, bmName)
								bmOutputDir = './{}/{}/'.format(gsOutDir, outputName)
								ElfinUtils.mkdir(bmOutputDir)

								configJson = OrderedDict([
									('inputFile', './{}/{}.json'.format(bmDir, bmName)),
									('xDBFile', './res/xDB.json'),
									('outputDir', bmOutputDir),
									('randSeed', '0x600d1337'),

									('chromoLenDev', cld),
									('gaPopSize', gps),
									('gaIters', gi),
									('gaSurviveRate', gsr),
									('gaCrossRate', gcr),
									('gaPointMutateRate', gpmr),
									('gaLimbMutateRate', glmr),
								])

								json.dump(configJson,
									open(gsConfDir + '/{}.json'.format(outputName), 'w'),
									separators=(',', ':'),
									ensure_ascii=False,
									indent=4)

							configId = configId + 1

	print 'Max GS config ID: {}'.format(configId - 1)

if __name__ == '__main__':
	main()
