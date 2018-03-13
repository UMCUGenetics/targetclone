class SimulationDataHandler:
	
	cMatrix = None
	aMatrix = None
	samples = None #This holds the chromosomes and positions
	afMatrix = None
	snvMatrix = None
	realTree = None
	chromosomes = None
	positions = None
	
	#After we finish running the simulations, we can use pickle to store all the simulation results in a file by uique identifier. This can be loaded by the permutation rounds
	
	def setSimulationData(self, cMatrix, aMatrix, samples, afMatrix, snvMatrix, realTree, chromosomes, positions):
		self.cMatrix = cMatrix
		self.aMatrix = aMatrix
		self.samples = samples
		self.afMatrix = afMatrix
		self.snvMatrix = snvMatrix
		self.realTree = realTree
		self.chromosomes = chromosomes
		self.positions = positions
		
		
		