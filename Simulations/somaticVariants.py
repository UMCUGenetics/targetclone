class SomaticVariant:
	
	chromosome = 0
	position = 0
	value = 0
	alleles = None
	ind = 0
	lost = False #keep track of if the variant has been lost before already, if true then we can never re-gain it due to the ISA. 
	
	def __init__(self, chromosome, position, value):
		self.chromosome = chromosome
		self.position = position
		self.value = value
		self.alleles = []
		
	def getChromosome(self):
		return self.chromosome
	
	def getPosition(self):
		return self.position
	
	def getValue(self):
		return self.value
	
	def getAlleles(self):
		return self.alleles