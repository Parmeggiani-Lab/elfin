#!/usr/bin/env python

import utils
import json
import argparse
import glob
import numpy as np
import GridSearchParams

# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt

def main():
	ap = argparse.ArgumentParser(description='Analyse Grid Search output files');
	# ap.add_argument('--gsConfDir', default='./gsConfigs/')
	# ap.add_argument('--gsOutDir', default='./gs_out/')
	ap.add_argument('--gsVersion', type=int, default=2)

	args = ap.parse_args()
	globals().update(vars(args))

	gsConfDir = './gsConfigsV{}/'.format(gsVersion)
	gsOutDir = './gsOutV{}/'.format(gsVersion)
	
	gsParams = utils.Bunch(GridSearchParams.getGSParams(gsVersion))
	
	scoreLineToken = '#0 score '
	aggregate = []
	confIdRange = xrange(0, gsParams.nRuns/len(gsParams.bmNames))
	for configId in confIdRange:
		for bmName in gsParams.bmNames:
			gsConfName = 'gs_{}_{}'.format(configId, bmName)
			gsConfFile = utils.normPath('{}/{}.json'.format(gsConfDir, gsConfName))
			gsOutLogFile = utils.normPath('{}/{}/log'.format(gsOutDir, gsConfName))

			conf = utils.readJSON(gsConfFile)

			with open(gsOutLogFile, 'r') as file:
				scoreLine = [l for l in file.read().split('\n') if scoreLineToken in l]

			tokenClearLine = scoreLine[0][scoreLine[0].find(scoreLineToken)+len(scoreLineToken):]
			colonClearLine = tokenClearLine[:tokenClearLine.rfind(':')]
			score = float(colonClearLine)
			aggregate.append({
				'gsConfName': gsConfName,
				'conf': conf,
				'score': score
				})

	# Average the 3 bm scores
	averaged = []
	for configId in confIdRange:
		gsConfHeader = 'gs_{}_'.format(configId)
		gsGroup = [a for a in aggregate if gsConfHeader in a['gsConfName']]

		assert(len(gsGroup) == 3)
		avgScore = sum([a['score'] for a in gsGroup]) / 3.0
		averaged.append({
			'gsConfName': gsGroup[0]['gsConfName'],
			'conf': gsGroup[0]['conf'],
			'score': avgScore
			})

	scores = [a['score'] for a in averaged]
	chromoLenDevs = [a['conf']['chromoLenDev'] for a in averaged]
	gaSurviveRates = [a['conf']['gaSurviveRate'] for a in averaged]
	gaCrossRates = [a['conf']['gaCrossRate'] for a in averaged]

	estimateIndependentBest(confIdRange,
							scores,
							averaged,
							chromoLenDevs,
							gaSurviveRates,
							gaCrossRates,
							gsParams.chromoLenDevs,
							gsParams.gaSurviveRates,
							gsParams.gaCrossRates,
							gsParams.pmRatios)

	estimateBestMedoid(confIdRange,
							scores,
							averaged,
							chromoLenDevs,
							gaSurviveRates,
							gaCrossRates,
							gsParams.chromoLenDevs,
							gsParams.gaSurviveRates,
							gsParams.gaCrossRates,
							gsParams.pmRatios)

	utils.pauseCode()

def estimateBestMedoid(confIdRange,
					   scores,
					   averaged,
					   chromoLenDevs,
					   gaSurviveRates,
					   gaCrossRates,
					   chromoLenDevsTicks,
					   gaSurviveRatesTicks,
					   gaCrossRatesTicks,
					   pmRatios):
 	scoreCutoff = 10.0
	eliteConfs = [a['conf'] for a in averaged if a['score'] < scoreCutoff]

	assert(len(eliteConfs) > 0)

	confKeys = ['chromoLenDev',
				'gaSurviveRate',
				'gaCrossRate',
				'gaPointMutateRate',
				'gaLimbMutateRate']
	rawMatrix = np.asarray([[el[k] for k in confKeys] for el in eliteConfs])

	# Normalise influences
	matMax = np.max(rawMatrix, 0)
	matMin = np.min(rawMatrix, 0)
	matRange = matMax - matMin
	for i in xrange(0, len(rawMatrix)):
		rawMatrix[i] = (rawMatrix[i] - matMin) / matRange;

	# Compute pairwise distances
	sumPwDists = []
	for r1 in rawMatrix:
		rowPwdSum = 0
		for r2 in rawMatrix:
			rowPwdSum += np.linalg.norm(r2 - r1)

		sumPwDists.append(rowPwdSum)

	chosenId = sumPwDists.index(min(sumPwDists))

	print 'Best medoid: {}'.format(json.dumps(eliteConfs[chosenId], sort_keys=True, indent=4))

def estimateIndependentBest(confIdRange,
							scores,
							averaged,
							chromoLenDevs,
							gaSurviveRates,
							gaCrossRates,
							chromoLenDevsTicks,
							gaSurviveRatesTicks,
							gaCrossRatesTicks,
							pmRatios):
	bestChromoLenDev 		 	= indepAvg(confIdRange, chromoLenDevsTicks, chromoLenDevs, scores)
	bestGASurviveRate 		 	= indepAvg(confIdRange, gaSurviveRatesTicks, gaSurviveRates, scores)
	bestGACrossRate 			= indepAvg(confIdRange, gaCrossRatesTicks, gaCrossRates, scores)

	# PM and LM rates depend on Cross rates
	bestGAPointMutateRates = []
	bestGALimbMutateRates = []
	gaCrossRateScores = []
	for cr in gaCrossRatesTicks:
		# Each remaining portion after cross rates
		# generate 3 ratios of RM and LM
		rem = 0.9999 - cr
		gaPointMutateRatesTicks 	= []
		gaLimbMutateRatesTicks  	= []
		for pmRatio in pmRatios:
			(pm, lm) = (pmRatio * rem, (0.9999 - pmRatio) * rem)
			gaPointMutateRatesTicks.append(pm)
			gaLimbMutateRatesTicks.append(lm)

		averagedCrGroup = [a for a in averaged if a['conf']['gaCrossRate'] == cr]
		gaPointMutateRates = [a['conf']['gaPointMutateRate'] for a in averagedCrGroup]
		gaLimbMutateRates = [a['conf']['gaLimbMutateRate'] for a in averagedCrGroup]

		idRange = xrange(0, len(averagedCrGroup))
		gaCrossRateScores.append(np.mean([acr['score'] for acr in averagedCrGroup]))
		bestGAPointMutateRates.append(indepAvg(idRange, gaPointMutateRatesTicks, gaPointMutateRates, scores))
		bestGALimbMutateRates.append(indepAvg(idRange, gaLimbMutateRatesTicks, gaLimbMutateRates, scores))

	print 'Indepdent variable average estimates:'
	print ('\tbest chromoLenDev={}\n' + \
		   '\tbest gaSurviveRate={}\n' + \
		   '\tbest gaCrossRates={}\n' + \
		   '\n\tgaCrossRatesTicks={}\n' + \
		   '\tgaCrossRateScores={}\n' + \
		   '\tbest gaPointMutateRates={}\n' + \
		   '\tbest gaLimbMutateRates={}\n').format(bestChromoLenDev,
		   											bestGASurviveRate,
		   											bestGACrossRate,
		   											gaCrossRatesTicks,
		   											gaCrossRateScores,
		   											bestGAPointMutateRates,
		   											bestGALimbMutateRates)

def indepAvg(confIdRange, ticks, confVals, scores):
	res = []
	for t in ticks:
		ids = [i for i in confIdRange if confVals[i] == t]

		filtered = [scores[i] for i in ids]
		assert(len(filtered) > 0)

		res.append(np.mean([scores[i] for i in ids]))
	return ticks[res.index(min(res))]

if __name__ == '__main__':
	utils.safeExec(main)