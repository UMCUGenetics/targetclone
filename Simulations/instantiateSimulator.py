import pickle

from simulations_permutations import Simulator


#Make TargetClone instance (later allow us to load the right instance based on kmin/kmax)
simulator = Simulator() #pre-initialize the method's C and mu combinations, required only once
simulator.initialize()

output = open('simulator.pkl', 'wb')

# Pickle dictionary using protocol 0.
pickle.dump(simulator, output)

output.close()