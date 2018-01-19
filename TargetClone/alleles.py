from editDistance import EditDistance

#Class to describe allele objects. Alleles know their A/B ratio and can output themselves in string format. 
class Alleles:
	ACount = None
	BCount = None
	alleleString = None #string representation of alleles
	alleleIdentifiers = None #unique identifier for the alleles, each allele has its own unique label (for simulations)

	def __init__(self, ACount, BCount):
		self.ACount = ACount
		self.BCount = BCount
		self.setAllelesAsString()
		self.alleleIdentifiers = []
	
	def getAllelesAsString(self):
		return ('A' * self.ACount) + ('B' * self.BCount)
	
	def setAllelesAsString(self): #re-run this every time the values of ACount and BCount update!
		self.alleleString = ('A' * self.ACount) + ('B' * self.BCount)
	
	#Compute the allele event distance between self and another set of alleles. 
	def computeEditDistance(self, otherAlleles):
		return EditDistance().editDistance(self.alleleString, otherAlleles.alleleString)