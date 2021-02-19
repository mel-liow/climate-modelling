##
# Ice sheet simulation
# A real-time evolving simulation of a 1-D ice sheet that models an 
# increasing profile elevation as snow accumulates over time. 
# This drives a pressure gradient which causes an outward flow in both directions.
##

import numpy as np
import matplotlib.pyplot as plt

# Predefined constants
nX = 10			# number of grid points
domainWidth = 1e6		# meters
timeStep = 100			# years
nYears = 20000			# years
flowParam = 1e4			# m horizontal / yr
snowFall = 1			# m / y
plotLimit = 4000

nSteps = int(nYears / timeStep)
dX = domainWidth / nX

# Initialise elevation and flow
elevations = np.zeros(nX+2)
flows = np.zeros(nX+1)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

# Loop through time
for i in range(0, nYears, timeStep):
	for ix in range(0, nX+1):
		surface_gradient = ( elevations[ix] - elevations[ix+1] ) / dX
		flows[ix] = surface_gradient * flowParam  * ( elevations[ix]+elevations[ix+1] ) / 2 / dX

	for ix in range(1, nX+1):
		elevations[ix]  = elevations[ix] + ( snowFall + flows[ix-1] - flows[ix] ) * timeStep


	print ("Years:", i)
	ax.clear()			# to update graph
	ax.plot(elevations)
	plt.title('1D Ice Sheet Model')
	plt.ylim(0,plotLimit)
	plt.show(block=False)
	plt.pause(0.001)			# delay between update

ax.clear()
ax.plot( elevations )
ax.set_ylim([0,plotLimit])
plt.show()