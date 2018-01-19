import numpy as np

#formally, this could inherit from a 'measurement' class. 
class LAF:

	measurements = None
	chromosomes = None
	starts = None
	ends = None
	segmentation = None
	
	def __init__(self, measurements, chromosomes, starts, ends):
		self.correctAfMeasurements(measurements)
		self.chromosomes = chromosomes
		self.starts = starts
		self.ends = ends
		self.cMuCombinations = []

	#Accept AF measurements, return LAF measurements
	def correctAfMeasurements(self, measurements):
		self.measurements = []
		for measurement in measurements:
			if measurement > 0.5:
				self.measurements.append(1-measurement)
			else:
				self.measurements.append(measurement)
				
	#THe segmentation applies to where we segment the measurements.
	def setSegmentation(self, segmentation):
		self.segmentation = segmentation
		
			
	