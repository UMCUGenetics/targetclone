import numpy as np

class Mu:
	#The given muIndex is always the mu value of the normal cell
	def __init__(self, muIndex, mu=None):
		if mu is None:
			self.mu=[]
			self.mu.append(muIndex/float(100))
			self.mu.append((100-muIndex)/float(100))
		else:
			self.mu = mu