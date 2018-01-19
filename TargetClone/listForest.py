#Process an input list of annotations for every iteration and plot the trees separately in tabs
import forestry
from bokeh.io import show,output_file,save

#The main class calls the function here that then using the information creates a top-10 of the trees. The tree visualization package is then invoked to generate output. 

def generateOutput(names, parents, annotations, weights, mu, messages, outputDir):

	#Do a sorting of the keys by ISA resolved & weights.
	#First obtain all the ISA resolved iterations, separate these from the ISA non-resolved ones
	#Sort both of the lists
	#Make a top-10 to report
	
	resolvedNames = dict()
	resolvedParents = dict()
	resolvedAnnotations = dict()
	resolvedWeights = dict()
	resolvedMu = dict()
	resolvedMessages = dict()
	
	unresolvedNames = dict()
	unresolvedParents = dict()
	unresolvedAnnotations = dict()
	unresolvedWeights = dict()
	unresolvedMu = dict()
	unresolvedMessages = dict()
	
	#Also get the total weight
	totalWeightResolved = []
	totalWeightUnresolved = []
	
	currentResolvedInd = 0
	currentUnresolvedInd = 0
	for iteration in range(0, len(names.keys())):
		if messages[iteration] == '':
			resolvedNames[currentResolvedInd] = names[iteration]
			resolvedParents[currentResolvedInd] = parents[iteration]
			resolvedAnnotations[currentResolvedInd] = annotations[iteration]
			resolvedWeights[currentResolvedInd] = weights[iteration]
			resolvedMu[currentResolvedInd] = mu[iteration]
			resolvedMessages[currentResolvedInd] = messages[iteration]
			
			totalWeightResolved.append(sum(weights[iteration]))
			currentResolvedInd += 1
		else:
			unresolvedNames[currentUnresolvedInd] = names[iteration]
			unresolvedParents[currentUnresolvedInd] = parents[iteration]
			unresolvedAnnotations[currentUnresolvedInd] = annotations[iteration]
			unresolvedWeights[currentUnresolvedInd] = weights[iteration]
			unresolvedMu[currentUnresolvedInd] = mu[iteration]
			unresolvedMessages[currentUnresolvedInd] = messages[iteration]
			
			totalWeightUnresolved.append(sum(weights[iteration]))
			currentUnresolvedInd += 1
		
	
	#Sort the resolved and unresolved iterations by weight
	sortedResolvedInd = [i[0] for i in sorted(enumerate(totalWeightResolved), key=lambda x:x[1])]
	sortedUnresolvedInd = [i[0] for i in sorted(enumerate(totalWeightUnresolved), key=lambda x:x[1])]


	#Combine both data into a top-10 list
	reportingTrees = 10
	if reportingTrees > len(names.keys()): #Check if we have enough to report, otherwise report fewer results
		reportingTrees = len(names.keys())
		
	#Make a combined set to get the trees from 
	top10Names = dict()
	top10Parents = dict()
	top10Annotations = dict()
	top10Weights = dict()
	top10Mu = dict()
	top10Messages = dict()
	
	currentResolvedInd = 0
	currentUnresolvedInd = 0
	for ind in range(0, reportingTrees):
		
		if ind >= len(sortedResolvedInd): #if we have more than the resolved ISAs, we start adding unresolved ones. 
			top10Names[ind] = (unresolvedNames[sortedUnresolvedInd[currentUnresolvedInd]])
			top10Parents[ind] = (unresolvedParents[sortedUnresolvedInd[currentUnresolvedInd]])
			top10Weights[ind] = (unresolvedWeights[sortedUnresolvedInd[currentUnresolvedInd]])
			top10Annotations[ind] = (unresolvedAnnotations[sortedUnresolvedInd[currentUnresolvedInd]])
			top10Mu[ind] = (unresolvedMu[sortedUnresolvedInd[currentUnresolvedInd]])
			top10Messages[ind] = (unresolvedMessages[sortedUnresolvedInd[sortedUnresolvedInd[currentUnresolvedInd]]])
			
			currentUnresolvedInd += 1
		else:
			top10Names[ind] = (resolvedNames[sortedResolvedInd[currentResolvedInd]])
			top10Parents[ind] = (resolvedParents[sortedResolvedInd[currentResolvedInd]])
			top10Weights[ind] = (resolvedWeights[sortedResolvedInd[currentResolvedInd]])
			top10Annotations[ind] = (resolvedAnnotations[sortedResolvedInd[currentResolvedInd]])
			top10Mu[ind] = (resolvedMu[sortedResolvedInd[currentResolvedInd]])
			top10Messages[ind] = (resolvedMessages[sortedResolvedInd[currentResolvedInd]])
			
			currentResolvedInd += 1
		
	#Show the trees in the final output

	#Loop through the lists and concatenate into different iteration lists
	allIterationData = []
	allIterationNames = []
	for iteration in range(0, len(top10Names.keys())):
		iterationNames = top10Names[iteration]
		iterationParents = top10Parents[iteration]
		iterationWeights = top10Weights[iteration]
		totalWeight = sum(iterationWeights)
		iterationAnnotations = top10Annotations[iteration]
		iterationMu = top10Mu[iteration]
		iterationMessages = top10Messages[iteration]
		
		
		iterationSet = [iterationNames,iterationParents,iterationMu, iterationWeights, iterationAnnotations, iterationMessages]
		allIterationData.append(iterationSet)
		allIterationNames.append(str(iteration+1) + ': ' + str(int(totalWeight))) #concatenate the names with the weights
		
	# optional:
	output_file(outputDir + '/trees.html', title='Inferred trees', mode='cdn', root_dir=None)
	# required:
	save(forestry.growForest(allIterationData,allIterationNames))

