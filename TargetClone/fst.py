import numpy as np
from editDistance import EditDistance
import re
import settings

#FST for computing event distances 
class FST:
	
	#COmpute distance based on copy numbers and report the distance. Profile 1 contains copy numbers of 1 subclone, profile 2 of another. This in general works between longer profiles or profiles with only one measurement position.  
	def computeDistance(self, profile1, profile2):
		#Given profile 1 and 2, go through the positions individually and increase the penalty for impossible situations
		
		totalDistance = 0
		previousAction = "match"
		
		for pInd in range(0, len(profile1)):
			
			if np.isnan(float(profile1[pInd])) or np.isnan(float(profile2[pInd])): #in this case we cannot compute a distance
				continue
				
			#Determine the action
			#The action is an amplification if sample - parent is positive
			#The action is a deletion if sample - parent is negative
			#The action is a match if sample - parent is 0
			action = "match"
			if (int(profile1[pInd]) - int(profile2[pInd]) > 0):
				action = "amp"
			elif (int(profile1[pInd]) - int(profile2[pInd]) < 0):
				action = "del"

			#if the action is not the same as the previous action, assign a penalty equal to the distance
			if action is not previousAction and action is not "match":
				totalDistance += abs(int(profile1[pInd]) - int(profile2[pInd])) #extra high penalty for larger distances. 
			#If we want, we can use the weight to determine if amplifications or deletions also have the same distance.
			
			previousAction = action
		return(totalDistance)
	
	#Computes P(C|T) based on the alleles of the subclones. THe parent alleles (pAlleles) and subclone alleles are provided for 2 positions independently. 
	def computePcBasedOnAlleles(self, pAlleles1, pAlleles2, alleles1, alleles2, eventDistance):
		#If the alleles are not the same as in the previous position,
		#Compute the event distance between the alleles (in total)
		#This is the combination score that we have. It should penalize preferring large differences in copy numbers for LOH regions.

		totalDistance = 0

		distancePenalty = settings.fst['pCWeight']
		#Check if the first alleles are different from the parent. The event distance between these is the distance.
		if alleles1 != pAlleles1:
			totalDistance += eventDistance.getRawEventDistance(pAlleles1, alleles1)
			

		#If then the second alleles are the same as the first, we return the previous distance (horizontal dependency).
		if alleles1 == alleles2:

			return totalDistance * distancePenalty
		else:
			if alleles2 != pAlleles2: #if the alleles are not the same, we have a new event and thus count new event distance. 
				totalDistance += eventDistance.getRawEventDistance(pAlleles2, alleles2)


		return totalDistance * distancePenalty #the penalty is higher horizontally, this helps in obtaining more equal copy numbers. 
	
	#Return all unique values in the input vectors. 
	def uniq(self, input):
		output = []
		for x in input:
			if x not in output:
				output.append(x)
		return output
	
	#This method is a gigantic mess, we should properly split it into neater functions and make it more usable
	#We check if we can find positions where we confidently see LOH in the parent but not in the child. If this is the case, then the parent cannot be the parent of the child. We check the same for alternating
	#alleles based on the allele frequency. 
	def computeAlleleDistance(self, profile1, profile2, sample1, sample2): #we need the LAF here of sample 1 (to determine if the LAF meets the LOH requirements), thus we provide the sample objects as well
		#Given profile 1 and 2, go through the positions individually and increase the penalty for impossible situations
		laf1 = sample1.measurements
		laf2 = sample2.measurements
			
		#Obtain all the settings that we use for the penalties
		previousAction = "match"
		lohLowerBound = settings.fst['lohPosLowerBound'] -1 
		notLohLowerBound = settings.fst['notLohPosLowerBound'] -1
		lafUpperBound = settings.fst['lafLohUpperBound']
		
		penalty = settings.fst['lossOfMediumConfidentLohPenalty']
		
		#Set the bounds for when we see one copy of one allele or one copy of the other based on the AF. 
		afLowerBound = lafUpperBound
		afUpperBound = 1 - lafUpperBound
		afPosLowerBound = settings.fst['alternatingAFPosLowerBound'] -1 #we need at least x adjacent positions to be LOH for the alternating allele checks to be used. 
		
		#if the previous position was LOH and we see at least x positions that are LOH, we will assign the LOH penalty.
		previousLoh = False
		#At each position, we need to check if the previous positions has a LAF smaller than the setting. If we have n consecutive positions matching this, it is really LOH. 
		#keep which positions are loh in the one sample. If these are not LOH in the other sample, then the distance is inf!
		currentConsecutiveAlternatingAf = 0
		currentConsecutiveLoh1 = 0 
		consecutiveLohRegion = []
		previousLoh1 = False
		currentConsecutiveLoh2 = 0
		lohFound1 = False
		lohFound2 = False
		alternatingAfFound = False
		messages = []
		previousChromosome = 1
		alleleALost = False
		alleleBLost = False
		alleleALostPos = 0
		alleleBLostPos = 0
		totalDistance = 0
		
		for pInd in range(0, len(profile1)):
			
			#if we have switched to the next chromosome, we reset the memory of the previous position.
			if previousChromosome != laf1.chromosomes[pInd]:
				previousAction = "match" #if we switch chromosome borders, we can never end up with a change where the previous action is the same.
				#reset the LOH counters when we cross chromosome borders
				lohFound1 = False
				lohFound2 = False
				currentConsecutiveLoh1 = 0
				currentConsecutiveLoh2 = 0
				parentNotLoh = 0
				childNotLoh = 0
				parentNotLohFound = False
				childNotLohFound = False
				alternatingAfFound = False
				currentConsecutiveAlternatingAf = 0
				alleleALost = False
				alleleBLost = False
				alleleALostPos = 0
				alleleBLostPos = 0

			#If a measurement is nan, we cannot compute a distance, so skip it. 
			if profile1[pInd] is np.nan or profile2[pInd] is np.nan:
				previousChromosome = laf1.chromosomes[pInd]
				continue
			if profile1[pInd] is None or profile2[pInd] is None:
				previousChromosome = laf1.chromosomes[pInd]
				continue
			
			#Check if we make a change or if it is a match. In case of a match, we are not interested in assigning penalties
			if profile1[pInd].ACount == profile2[pInd].ACount and profile1[pInd].BCount == profile2[pInd].BCount:
				action = "match"
			else:
				action = "change"

			#set if we find LOH at this position in the parent or child in our A matrix predictions.
			parentLoh = False
			childLoh = False
			if profile1[pInd].ACount < 1 and profile1[pInd].BCount > 0:
				parentLoh = True
			if profile2[pInd].ACount < 1 and profile2[pInd].BCount > 0:
				childLoh = True
			
			#If both the measurements and the method report LOH
			methodMeasurementsLohReportedParent = False
			methodMeasurementsLohReportedChild = False
			
			#Here we check if we confidently find LOH in the parent (both in the measurements and method), or only in the measurements, or not at all. We reset the counters for if we see LOH based on these choices
			if laf1.measurements[pInd] <= lafUpperBound and parentLoh is True and previousChromosome == laf1.chromosomes[pInd]:
				currentConsecutiveLoh1 += 1
				
				#if we see in the measurements of the parent that there is clearly LOH, we reset the not LOH counters. 
				parentNotLoh = 0 #we ignore regions where only the method reports LOH, this is not confident enough for annotations. 
				parentNotLohFound = False
				methodMeasurementsLohReportedParent = True
			elif laf1.measurements[pInd] <= lafUpperBound and parentLoh is False and previousChromosome == laf1.chromosomes[pInd]:
				#we end up here when the measurements are low, but the method thinks the region is not LOH.
				#we consider these regions as a penalty. 
				if action is not previousAction and action is not "match":
					totalDistance += penalty
			else: #if we end up here, we are in a region where the parent is not LOH. 
				#save if we introduced new LOH only if we are confident.
				
				#If we end up back at a position where the parent does not have LOH and we have not found 10 LOH positions up until now in the parent,
				#this is a new event.
				
				if lohFound2 is True and parentNotLohFound is True:
					affectedRegion = laf1.segmentation.getSegmentName(laf1.chromosomes[pInd], laf1.starts[pInd], laf1.ends[pInd])
					consecutiveLohRegion.append(affectedRegion)
					
				currentConsecutiveLoh1 = 0
				lohFound1 = False
				parentNotLoh += 1
				#consecutiveLohRegion = []
			
			#Repeat LOH detection for the child
			#There is only LOH in the child if the measurements support it and the method thinks that as well (avoid recognizing ie ABB as LOH when there is low normal cell contamination). 
			if laf2.measurements[pInd] <= lafUpperBound and childLoh is True and previousChromosome == laf2.chromosomes[pInd]:
				currentConsecutiveLoh2 += 1
				childNotLoh = 0
				childNotLohFound = False
				methodMeasurementsLohReportedChild = True
				
			elif laf2.measurements[pInd] > lafUpperBound and previousChromosome == laf2.chromosomes[pInd]: #if the measurements do not support LOH, but the method does, there is no LOH. 
				currentConsecutiveLoh2 = 0
				lohFound2 = False
				childNotLoh += 1 #we are not confident enough based on the method alone that this region is LOH.
				
			else: #we only end up here if the method does not report LOH and the measurements are also high. 
				#if we find confident LOH in the child that is not in the parent, it is a new event.
				#These events are now not correctly added 

				currentConsecutiveLoh2 = 0
				lohFound2 = False
				childNotLoh += 1
			
			#check all positions, if we find 10 consecutive loh positions in sample 1, then this is loh. we check the same positions in sample 2. if we reach the
			#end of the loh in sample 1 and we do not find 10 consecutive in sample 2, this relation is not possible. 

			if currentConsecutiveLoh1 >= lohLowerBound:
				lohFound1 = True
				
			if currentConsecutiveLoh2 >= lohLowerBound:
				lohFound2 = True
				
			if parentNotLoh >= notLohLowerBound:
				parentNotLohFound = True
			if childNotLoh >= notLohLowerBound:
				childNotLohFound = True
			
			#as soon as we are sure that we have found confident loh, restrict the relations. 
			if lohFound1 is True and childNotLohFound is True:
		
				return [messages, settings.fst['lossOfConfidentLohPenalty']]
			
			#we only increase our alternating alleles found when we see LOH in both the child and parent at this position. 
			if methodMeasurementsLohReportedChild is True and methodMeasurementsLohReportedParent is True:
				#check if the AFs alternate
				if (sample1.afMeasurements[pInd] < afLowerBound and sample2.afMeasurements[pInd] > afUpperBound) or (sample2.afMeasurements[pInd] < afLowerBound and sample1.afMeasurements[pInd] > afUpperBound):
					currentConsecutiveAlternatingAf += 1
					
				
			#If we have enough evidence for alternating alleles, restrict the relations between the parent and child. 
			if currentConsecutiveAlternatingAf >= afPosLowerBound:
				alternatingAfFound = True
				#if we find alternating AF for at least 4 positions, we restrict the relationship between the samples.
			
				return [messages, settings.fst['alternatingAFPenalty']]
			
			#This is not used, but could be used to annotate in the tree if allele A or B has been lost to show if a different parental allele is lost in specific samples. 
			if methodMeasurementsLohReportedChild is True: 
				#Check if we should annotate if allele A or B has been lost
				if sample2.afMeasurements[pInd] < afLowerBound:
					alleleALostPos += 1
				if sample2.afMeasurements[pInd] > afUpperBound:
					alleleBLostPos += 1
				
				if alleleALostPos > lohLowerBound:
					alleleALost = True
				if alleleBLostPos > lohLowerBound:
					alleleBLost = True
			
			#if the parent has LOH in both the measurements and the child does not, restrict the region.
			#if the parent has LOH in the measurements and the child does not, penalize. 

			#if the action is not the same as the previous action, assign a penalty equal to the event distance

			if action is not previousAction and action is not "match":
				
				totalDistance += (EditDistance().editDistance(profile1[pInd].getAllelesAsString(), profile2[pInd].getAllelesAsString())) #It may be useful to provide the pre-computed distances here. Unless it's fast enough. 
		
			previousAction = action
			previousChromosome = laf1.chromosomes[pInd]
		
		#Set the events that we will annotate in the tree. 
		#if the set of LOH regions contains both the p and q arm, we only add the chromosome name.
		if len(self.uniq(consecutiveLohRegion)) > 0:
			#all chromosomes are added at once, so we should separate these
			#if we find two chromosome arms for the same chromosome, we merge these.
			#split each unique element after the number. If there are more than 2 of some number, this is the whole chromosome.
			#otherwise, we add the arm itself.
			chromosomes = []
			for element in self.uniq(consecutiveLohRegion):
				chromosomes.append(re.compile("\d+").search(element).group())
			
			#if an element has already been seen, we add the whole chromosome number.
			#then, we look at the unique set again and add the ones that have not been added.
			#for every chromosome number,
			#check if the consecutive loh regions contain that number with a "p" or a "q".
			mergedChromosomeList = []
			for chromosome in self.uniq(chromosomes):
				pFound = False
				qFound = False
				pPattern = re.compile(chromosome + "p")
				qPattern = re.compile(chromosome + "q")
				for element in self.uniq(consecutiveLohRegion):
					if pPattern.search(element):
						pFound = True
					if qPattern.search(element):
						qFound = True
				if pFound is True and qFound is True:
					mergedChromosomeList.append(chromosome)
					
					messages.append('+ LOH ' + chromosome)
				elif pFound is True and qFound is False:
					mergedChromosomeList.append(chromosome + "p")
					
					messages.append('+ LOH ' + chromosome + "p")
				elif pFound is False and qFound is True:
					mergedChromosomeList.append(chromosome + "q")
					
					messages.append('+ LOH ' + chromosome + "q")
		return [messages, totalDistance]
		
	