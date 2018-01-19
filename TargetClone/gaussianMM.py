import numpy as np
from scipy.stats import norm
from scipy.signal import argrelextrema
import math

#Class to generate our gaussian mixture model to obtain the probability that a LAF is measured given a C mu combination. 
class GaussianMM:

	mu = None
	sigma = None
	ratios = None
	pdf = None #Will hold the pdf
		
	def __init__(self, mu, sigma, ratios):
		self.mu = mu
		self.sigma = sigma
		self.ratios = ratios
		self.computePdf()
	
	#Compute the mixture model probabilities
	def computePdf(self):
		self.pdf = 0

		#Obtain all the mu values of the distribution, these are the LAF values that we can measure with a C mu combination
		uniqueMu = np.unique(np.asarray(self.mu))
		npRatios = np.asarray(self.ratios) #The ratios are the heights of the components, these are based on the probabilities
		uniqueRatios = []
		for i in range(0, len(uniqueMu)):
			indexes = [j for j,x in enumerate(self.mu) if x == uniqueMu[i]]
			uniqueRatios.append(np.sum(npRatios[indexes]))
		#Compute the density. First for the first component, then for every next LAF we again make a component and we multiply all together to make a Gaussian MM. 
		density = self.ratios[len(self.ratios)-1] * norm(self.mu[len(self.mu)-1], self.sigma).pdf(np.arange(0,1,0.01))
		for i in range(0, len(self.mu)-1):
			density += self.ratios[i] * norm(self.mu[i], self.sigma).pdf(np.arange(0,1,0.01))
		
		#Normalize the density between 0 and 1
		x = 0
		y = max(uniqueRatios)
		minDensity = min(density)
		maxDensity = max(density)
		densityRange = maxDensity - minDensity
		normalizedDensity = (density - minDensity) / float(densityRange)
		
		densityRatioRange = y - x
		normalizedDensity = (normalizedDensity*densityRatioRange) + x
		
		#Trim the density at 0.5
		self.pdf = normalizedDensity
		self.pdf[range(51,100)] = 0
		
	#Obtain a probability from the pdf at position x	
	def getPdf(self, x):
		return self.pdf[int(x*100)]
	
	#Compute the threshold locations
	#For the MM, mu means LAF, not actual mu...
	def getValleys(self):
		#Valleys are not working properly, due to the noise we sometimes end up with having no proper valleys if the laf values are close together!
		
		#Make the mu unique, make sure that we use LAF rather than AF to avoid problems further down the code
		#It is better to not do this with numpy, that takes much longer
		uniqueMu = []
		for mu in self.mu:
			if mu not in uniqueMu:
				uniqueMu.append(mu)
		
		
		#What we can instead try is to find the laf we are closest to.
		#for every mu, compute the mean compared to the previous entry. Then the thresholds will be at the mean between every adjacent pair of LAF in the distribution. 
		muAverages = []
		for muInd in range(1, len(uniqueMu)):
		 	avg = (uniqueMu[muInd] + uniqueMu[muInd-1]) / float(2)
		 	muAverages.append(round(avg,2))
		
		return muAverages
		