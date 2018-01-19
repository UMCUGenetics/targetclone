import pickle

from run import TargetClone


#Make TargetClone instance (later allows us to load the right instance based on kmin/kmax, this saves a lot of time because we do not need to make all combinations of C and mu each time, these can be re-used)
targetClone = TargetClone() #pre-initialize the method's C and mu combinations, required only once

output = open('targetClone.pkl', 'wb')

# Pickle dictionary using protocol 0.
pickle.dump(targetClone, output)

output.close()