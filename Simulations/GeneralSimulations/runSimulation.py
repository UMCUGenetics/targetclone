#!/usr/bin/python

import time
from simulation_generic import Simulator
import sys

#THis script will be called from a bash script to run on the HPC.

#It should load an instantiated simulator, and use this information to run a simulation 100 times. (Maybe we will start with running the simulator only 50 times, depending on how long it takes to run 1
#round of 100 permutations on the HPC).

#We need to:
#1. Make sure that the simulator is instantiated and saved (using pickle) -> is this necessary? It should not be that slow. 
#2. Make sure that we can run 1 simulation round if we call the expand subclone function
#3. Make sure that the data is written to files correctly, where:
#	-> Every simulation is stored in its own folder with an unique identifier. Here we obtain the underlying used data (no permutations)
#	-> The error is written to a new file (or will appending also work?)
#4. After the simulations have run, we can combine the errors together and make the boxplots and confusion matrix plots.

taskId = int(sys.argv[1]) -1 #we use the task ID as mu value, 
uniqueID = sys.argv[2]
#leaveOut = sys.argv[3]

#2.
start_time = time.time()
simulator = Simulator()
simulator.initializeGenericSimulation(taskId, uniqueID)
#simulator.subclonalExpansion()
print("--- %s seconds for 1 simulation run with 10 permutations ---" % (time.time() - start_time))
#For running a simulation, we want to have a number of parameters that we can set beforehand. These are for example an increase in noise level, and also the mu range to select from.
#THe noise level can be provided to the script (initialize). The mu range can also be defined there. The script should always select a random mu from this range!
#we then provide the parameters to this script, such that we have 10 jobs where the mu is in a certain range, etc. For the noise levels we can du 