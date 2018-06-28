# This file contains versioned grid search parameter data
def getGSParams(gsVersion):
	if gsVersion == 1:
		return getGSParamsV1()

	if gsVersion == 2:
		return getGSParamsV2()

	print 'Unknown Grid Search version: {}'.format(gsVersion)

def getGSParamsV2():
	# Define the grid
	chromoLenDevs 		= [0.2]
	gaPopSizes 			= [int(100e3)] 					
	gaIters 			= [int(1e3)]						
	gaSurviveRates 		= [0.02, 0.03, 0.05]
	gaCrossRates 		= [0.2, 0.4, 0.6]
	gaPointMutateRates 	= []
	gaLimbMutateRates  	= []

	# Create 3 configs - for 3 different benchmarks shapes
	# using the same config
	bmNames = ['6vjex8d', '9y8hxgo', 'j0m06n4']

	# PM and LM rates depend on Cross rates
	pmRatios = (0.3, 0.5, 0.70)
	for cr in gaCrossRates:
		# Each remaining portion after cross rates
		# generate 3 ratios of RM and LM
		rem = 0.9999 - cr
		for pmRatio in pmRatios:
			(pm, lm) = (pmRatio * rem, (0.9999 - pmRatio) * rem)
			gaPointMutateRates.append(pm)
			gaLimbMutateRates.append(lm)

	nRuns = len(chromoLenDevs) * len(gaPopSizes) * len(gaIters) * \
		len(gaSurviveRates) * len(gaCrossRates) * len(gaPointMutateRates) * \
		len(bmNames)

	return {
		'chromoLenDevs': chromoLenDevs,
		'gaPopSizes': gaPopSizes,
		'gaIters': gaIters,
		'gaSurviveRates': gaSurviveRates,
		'gaCrossRates': gaCrossRates,
		'gaPointMutateRates': gaPointMutateRates,
		'gaLimbMutateRates': gaLimbMutateRates,
		'bmNames': bmNames,
		'pmRatios': pmRatios,
		'nRuns': nRuns
	}

def getGSParamsV1():
	# Define the grid
	chromoLenDevs 		= [0.1, 0.2, 0.3]
	gaPopSizes 			= [int(100e3)] 					
	gaIters 			= [int(1e3)]						
	gaSurviveRates 		= [0.005, 0.01, 0.02]
	gaCrossRates 		= [0.3, 0.5, 0.7]
	gaPointMutateRates 	= []
	gaLimbMutateRates  	= []

	# Create 3 configs - for 3 different benchmarks shapes
	# using the same config
	bmNames = ['6vjex8d', '9y8hxgo', 'j0m06n4']

	# PM and LM rates depend on Cross rates
	pmRatios = (0.25, 0.5, 0.75)
	for cr in gaCrossRates:
		# Each remaining portion after cross rates
		# generate 3 ratios of RM and LM
		rem = 0.9999 - cr
		for pmRatio in pmRatios:
			(pm, lm) = (pmRatio * rem, (0.9999 - pmRatio) * rem)
			gaPointMutateRates.append(pm)
			gaLimbMutateRates.append(lm)

	nRuns = len(chromoLenDevs) * len(gaPopSizes) * len(gaIters) * \
		len(gaSurviveRates) * len(gaCrossRates) * len(gaPointMutateRates) * \
		len(bmNames)

	return {
		'chromoLenDevs': chromoLenDevs,
		'gaPopSizes': gaPopSizes,
		'gaIters': gaIters,
		'gaSurviveRates': gaSurviveRates,
		'gaCrossRates': gaCrossRates,
		'gaPointMutateRates': gaPointMutateRates,
		'gaLimbMutateRates': gaLimbMutateRates,
		'bmNames': bmNames,
		'pmRatios': pmRatios,
		'nRuns': nRuns
	}