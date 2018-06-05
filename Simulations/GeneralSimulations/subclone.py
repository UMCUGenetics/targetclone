import sys
sys.path.insert(0, '../../TargetClone/')


import numpy as np
import re
import random
from c import C
from alleles import Alleles
from copy import deepcopy
import uuid
from random import randint


class Subclone:
	
	C = None
	A = None
	somaticVariants = None
	parent = None
	children = None #we make a binary expansion, so there will always be only 1 new child, the other branch is the clone itself. 
	name = None
	
	def __init__(self):
		self.children = []
	
	def expand(self, duplicate):
		1+1
		
	def duplicateGenome(self, simulator):
		newC = []
		for c in self.C:
			oldCTumor = c.c[1]
			newC.append(C([2, oldCTumor*2]))
		newAlleles = []
		#Make sure that the allele identifiers are the same across a region with the same chromosome arm.
		
		prevArm = simulator.chromosomeArms[0]
		currentAlleleA = 'A' + str(uuid.uuid4())
		currentAlleleB = 'B' + str(uuid.uuid4())
		for a in range(0, len(self.A)):
			newACount = self.A[a].ACount*2
			newBCount = self.A[a].BCount*2
			#make the new names
			
			#Check for each arm if it is the same as the previous. If true, append the same name
			#if the chromosome arm is different, append a new one
			if simulator.chromosomeArms[a] != prevArm:
				currentAlleleA = 'A' + str(uuid.uuid4())
				currentAlleleB = 'B' + str(uuid.uuid4())
			
			#the identifiers for the new alleles depends on the ACount.
			newAlleleObjects = Alleles(newACount, newBCount)
			newAlleleObjects.alleleIdentifiers += self.A[a].alleleIdentifiers
			for newAIdentifier in range(0, (self.A[a].ACount*2 - self.A[a].ACount)):
				newAlleleObjects.alleleIdentifiers.append(currentAlleleA)
			for newBIdentifier in range(0, (self.A[a].BCount*2 - self.A[a].BCount)):
				newAlleleObjects.alleleIdentifiers.append(currentAlleleB)	
			newAlleles.append(newAlleleObjects)
			prevArm = simulator.chromosomeArms[a]
		
		#Check if our alleles always match across the chromosomes
		
		#we are using this for duplicating a precursor that does not have somatic variants. If the precursor did have somatic variants, we should duplicate these too here!
		
		precursor = Subclone()
		precursor.C = newC
		precursor.A = newAlleles
		precursor.somaticVariants = deepcopy(self.somaticVariants)
		precursor.parent = self
		precursor.name = str(uuid.uuid4())
		self.children.append(precursor)
		
		return precursor

	def expandSubclone(self, simulationProbabilities, simulator): #simulation probabilities can be accessed through simulator, this is not necessary to provide!
		
		
	
		
		#Apply the gains and losses
		
		newSubclone = Subclone()
		newSubclone.C = deepcopy(self.C)
		newSubclone.A = deepcopy(self.A)
		newSubclone.somaticVariants = deepcopy(self.somaticVariants)
		newSubclone.parent = deepcopy(self)
		self.children.append(newSubclone)
		newSubclone.name = str(uuid.uuid4())
		
		
		##This part can go to the probability handler, we duplicate it now
		possibleArms = []
		armLossProbabilities = []
		for armProbabilityList in simulationProbabilities.armProbabilities:
			possibleArms.append(armProbabilityList[0])
			armLossProbabilities.append(float(armProbabilityList[2]))
		
		
		#convert probabilities to numpy array	
		armLossProbabilities = np.asarray(armLossProbabilities)
		
		#These probabilities are not normalized, do this too:
		normalizedArmLossProbabilities = (armLossProbabilities / sum(armLossProbabilities))
		
		filteredPossibleArms = []
		filteredNormalizedArmLossProbabilities = []
		for arm in range(0, len(possibleArms)):
			chrInd = simulator.chromosomeArms.index(possibleArms[arm])
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == possibleArms[arm]]
			if newSubclone.C[chrIndices[0]].c[1] > (simulator.kmin-1): #we allow the loss of an arm, but if the total number of copies drops below kmin, the subclone is not viable.  
				filteredPossibleArms.append(possibleArms[arm])
				filteredNormalizedArmLossProbabilities.append(armLossProbabilities[arm])
	
		#If we remove some arms from this list, we also need to make sure that we re-normalize the probabilities. 		
		filteredNormalizedArmLossProbabilities = filteredNormalizedArmLossProbabilities / sum(filteredNormalizedArmLossProbabilities)
		
		
		#Determine how many changes 
		maxChanges = len(possibleArms) #Every arm has a possibility of being hit by a gain or loss once.
		
		#1. Choose a random number of changes to make
		numChanges = randint(0, (maxChanges-1))
		
		#2. Randomly divide these into gains and losses
		division = randint(0, numChanges)
		
		armGainCount = numChanges - division
		armLossCount = division
		
		#Determine how many whole chromosomes will be affected, again separate in gains and losses (are we even using this?)
		maxWcChanges = 22
		numWcChanges = randint(0, (maxChanges-1))
		
		wcDivision = randint(0, numWcChanges)
		
		wcGainCount = numWcChanges - wcDivision
		wcLossCount = wcDivision
		
		#Determine how many SNVs are gained
		snvGainCount = randint(0, (simulator.snvNum-1))
		
		
		##Determine how many somatic variants we have per allele (the alleles should be the same across all indices of the chromosome arms)
		allelesAndVariantCounts = dict()
		for arm in possibleArms:
			armIndex = possibleArms.index(arm)
			chrInd = simulator.allChromosomeArms.index(arm)
			chrIndices = [i for i, x in enumerate(simulator.allChromosomeArms) if x == arm]
			for variant in newSubclone.somaticVariants:
				#Get the somatic variants that are present on this arm
				identifiers = variant.alleles
				if variant.chromosome == arm:
					#get the alleles
					
					for identifier in identifiers:
						if identifier not in allelesAndVariantCounts.keys():
							allelesAndVariantCounts[identifier] = []
							
						allelesAndVariantCounts[identifier].append(variant)
				else:
					
					for identifier in newSubclone.A[chrIndices[0]].alleleIdentifiers:
						allelesAndVariantCounts[identifier] = []
		
		probableElements = np.where(filteredNormalizedArmLossProbabilities > 0)[0]
		#print "pe: ", probableElements
		#Find how many elements in the normalized array are larger than 0. 
		armsToSample = 0
		if len(probableElements) < armLossCount:
			armsToSample = len(probableElements)
		else:
			armsToSample = armLossCount
			
		armsLost = np.random.choice(filteredPossibleArms, armsToSample, p=filteredNormalizedArmLossProbabilities, replace=False)

		#Update the lost arms in the C and A matrix
		for arm in range(0, len(armsLost)):
			#we need to find the right index, arms lost contains the values
			lostArmInd = possibleArms.index(armsLost[arm])
			
			
			chrInd = simulator.allChromosomeArms.index(possibleArms[lostArmInd])
			chrIndices = [i for i, x in enumerate(simulator.allChromosomeArms) if x == possibleArms[lostArmInd]]
			
			
			
			
			#From all alleles, check which one has the least somatic variants. If there is a tie, we choose a random one.
			leastSVAllele = ''
			currentBest = (simulator.snvNum + 1) #we can never has less than the maximum number of SVs
			#print "we have identifiers: ", newSubclone.A[lostArmInd].alleleIdentifiers
			#if there is no identifier, we have no alleles left and the clone is not viable.
			if len(newSubclone.A[chrIndices[0]].alleleIdentifiers) < 1:
				return False
			
			for allele in newSubclone.A[chrIndices[0]].alleleIdentifiers:
				#check how many somatic variants the alleles have
				if len(allelesAndVariantCounts[allele]) < currentBest:
					leastSVAllele = allele
				#randomAllele = random.sample(newSubclone.A[lostArmInd].alleleIdentifiers, 1)[0] #select a random allele to be lost
			
			#here we wish to remove the allele with the least number of somatic variants
				
			randomAlelleIndex = newSubclone.A[chrIndices[0]].alleleIdentifiers.index(leastSVAllele)
		
			for ind in chrIndices:
				
				newSubclone.C[ind].c[1] = newSubclone.C[ind].c[1] - 1 #only lose the copy after we are sure that we can lose the allele as well
				
				
				
				del newSubclone.A[ind].alleleIdentifiers[randomAlelleIndex]
				#if the random allele matches on A, we remove the ACount. Otherwise, reduce the Bcount.
	
				if re.match('^A.+', leastSVAllele) is not None:
					newSubclone.A[ind].ACount = newSubclone.A[ind].ACount - 1 #we always keep allele B
					#also remove the selected allele from the list
				if re.match('^B.+', leastSVAllele) is not None:
					newSubclone.A[ind].BCount = newSubclone.A[ind].BCount - 1 #we always keep allele A
			
			lostVariants = allelesAndVariantCounts[leastSVAllele]
			
			varCounter = 0
			for lostVariant in lostVariants:
				lostVariant.lost = True
				lostVariant.value = 0
				varCounter += 1
	
		#Gain alleles: it is always possible to gain, but only to gain alleles that have not yet been lost
		
		#The alleles to gain also depends on which alleles can still be gained.
		filteredPossibleArms = []
		filteredNormalizedArmLossProbabilities = []
		for arm in range(0, len(possibleArms)):
			chrInd = simulator.allChromosomeArms.index(possibleArms[arm])
			chrIndices = [i for i, x in enumerate(simulator.allChromosomeArms) if x == possibleArms[arm]]
			if newSubclone.C[chrIndices[0]].c[1] < (simulator.kmax+1): #we allow the loss of an arm, but if the total number of copies drops below kmin, the subclone is not viable.  
				filteredPossibleArms.append(possibleArms[arm])
				filteredNormalizedArmLossProbabilities.append(armLossProbabilities[arm])
		
		#If we remove some arms from this list, we also need to make sure that we re-normalize the probabilities. 		
		filteredNormalizedArmLossProbabilities = filteredNormalizedArmLossProbabilities / sum(filteredNormalizedArmLossProbabilities)
		
		#If we gain the allele that contains a somatic variant, we need to duplicate these as well. 
		armsGained = np.random.choice(filteredPossibleArms, armGainCount, p=filteredNormalizedArmLossProbabilities, replace=False)
		
		for arm in range(0, len(armsGained)):
			#The selected allele needs to be the same for all chromosomes. We only gain the allele once and always the same type!
			
			gainedArm = possibleArms.index(armsGained[arm])
			chrInd = simulator.allChromosomeArms.index(armsGained[arm])
			chrIndices = [i for i, x in enumerate(simulator.allChromosomeArms) if x == armsGained[arm]]
			
			
			aAlleles = []
			bAlleles = []

			for allele in newSubclone.A[chrIndices[0]].alleleIdentifiers:
				#if the arm matches on A, we add it to the A alleles, otherwise the B alleles.
				
				if re.match('^A.+', allele) is not None:
					aAlleles.append(allele)
					#print "a alleles: ", aAlleles
				else:
					bAlleles.append(allele)
					#print "b alleles: ", bAlleles
			
			#Here we do not allow an allele to be gained when the allele has already been lost!
			newAllelePrefix = 'A'
			if newSubclone.A[chrIndices[0]].ACount < 1 and newSubclone.A[chrIndices[0]].BCount > 0:
				for ind in chrIndices:
					newSubclone.A[ind].BCount = newSubclone.A[ind].BCount + 1
				#also choose which one to gain out of the B alleles
				randomAllele = random.sample(bAlleles, 1)[0]
				#then when we have selected one, we duplicate it in the allele identifiers
				newAllelePrefix = 'B'
			elif newSubclone.A[chrIndices[0]].BCount < 1 and newSubclone.A[chrIndices[0]].ACount > 0:
				for ind in chrIndices:
					newSubclone.A[ind].ACount = newSubclone.A[ind].ACount + 1
				randomAllele = random.sample(aAlleles, 1)[0]
			elif newSubclone.A[chrIndices[0]].BCount > 0 and newSubclone.A[chrIndices[0]].ACount > 0: #this should generally always be the case if we make sure that the kmin is not exceeded
				#select any old allele
				#why is there sometimes no allele here? Either of the above should work. This only means that sometimes even
				#print "identifier: ", newSubclone.A[gainedArm].alleleIdentifiers
				randomAllele = random.sample(newSubclone.A[chrIndices[0]].alleleIdentifiers, 1)[0]
				#print randomAllele
				if re.match('^A.+', randomAllele) is not None:
					for ind in chrIndices:
						newSubclone.A[ind].ACount = newSubclone.A[ind].ACount + 1 #we always keep allele B
				if re.match('^B.+', randomAllele) is not None:
					for ind in chrIndices:
						newSubclone.A[ind].BCount = newSubclone.A[ind].BCount + 1 #we always keep allele A
					newAllelePrefix = 'B'
			else: #if we have no copies left, we cannot do anything and w should skip this gain.
				continue
			newAlleleName = newAllelePrefix + str(uuid.uuid4())
			for ind in chrIndices:
				newSubclone.A[ind].alleleIdentifiers.append(newAlleleName)
			#we also notify the somatic variant that it has been duplicated
			duplicatedVariants = allelesAndVariantCounts[randomAllele]
			for duplicatedVariant in duplicatedVariants:
				duplicatedVariant.alleles.append(newAlleleName)
				
			#Only update C when we get here
			for ind in chrIndices:
			
				
				newSubclone.C[ind].c[1] = newSubclone.C[ind].c[1] + 1
			
		
		#Gain somatic variants: only gain somatic variants that have not been lost yet! We need to give these a tag to see if they have been lost or not in a previous iteration.
		#Currently losses are not implemented yet, so we can go with everything that is not already there.
		
		availableSomVarSlots = []
		for somVar in newSubclone.somaticVariants:
			#also check if the allele is still available, otherwise we cannot gain an allele there. 
			chromosome = somVar.chromosome
			chromosomeInd = possibleArms.index(somVar.chromosome)
			chrInd = simulator.allChromosomeArms.index(somVar.chromosome)
			chrIndices = [i for i, x in enumerate(simulator.allChromosomeArms) if x == somVar.chromosome]
			#print "the content: ", newSubclone.A[chromosomeInd].alleleIdentifiers
			if len(newSubclone.A[chrIndices[0]].alleleIdentifiers) < 1: #if there are no alleles to add the variant to, we skip it. 
				continue
			
			if somVar.value < 1 and somVar.ind not in simulator.usedVariantSlots: #Make sure that the position has not been used already.
				availableSomVarSlots.append(somVar)
		
		
		somaticVariantsGained = []
		if len(availableSomVarSlots) == 1: #if we cannot collect 2 variants to gain but at least one, we also need to sample only one
			somaticVariantsGained = random.sample(availableSomVarSlots, 1)
		if len(availableSomVarSlots) > 1:
			somaticVariantsGained = random.sample(availableSomVarSlots, snvGainCount)
		
		if len(somaticVariantsGained) > 0: #if we can gain anything at all
			#somaticVariantsGained = random.sample(availableSomVarSlots, simulator.cycleSVGains)
			
			for somaticVariantGained in somaticVariantsGained:
				somaticVariantGained.value = 1
				simulator.usedVariantSlots.append(somaticVariantGained.ind)
				#it knows on which chromosome it has been gained on, but we still need to assign an allele identifier
				#if the chromosome arm already has alleles, we can use these
				#we consider the entire arm to consist of one big allele
				variantArm = somaticVariantGained.chromosome
				variantArmInd = possibleArms.index(variantArm)
				chrIndices = [i for i, x in enumerate(simulator.allChromosomeArms) if x == variantArm]
				#Check the alleles that it has and select a random one
				randomAllele = random.sample(newSubclone.A[chrIndices[0]].alleleIdentifiers, 1)[0]
				somaticVariantGained.alleles.append(randomAllele)
		
		gainedPositions = []
		for var in somaticVariantsGained:
			gainedPositions.append(var.ind)
		
		
		
		#Also check for viability, report a message to the method when we cannot go any further from here because the subclone is not viable (a bit of a recursive fashion).
		#Our viability checks are:
		#- We may not exceed our kmin and kmax
		#- 12p amplifications are not allowed to be lost (is this possible now?)
		i = 0
		for cInd in range(0, len(newSubclone.C)):
			c = newSubclone.C[cInd]
		
			if c.c[1] < simulator.kmin or c.c[1] > simulator.kmax:
				return False
			i += 1
		
		# for a in range(0, len(newSubclone.A)):
		# 	print "chromosome arm: ", simulator.chromosomeArms[a]
		# 	print "alleles: ", newSubclone.A[a].alleleIdentifiers
		# exit()	
		return newSubclone
		