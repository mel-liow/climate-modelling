##
# Albedo Hysteresis loop
# An iterative model of the Albedo feedback and how it affects
# the temperature of the Earth.
##

import numpy as np
import matplotlib.pyplot as plt

# Predefined constants
epsilon = 1 									# Planet is assumed to be a black body
sigma = 5.67E-8								# Stafan-Boltzman constant (W/m^2K^4)
albedo_gradient = -0.01   		# Gradient of albedo-temperature relationship (K-1)
albedo_intercept = 2.8

# Predefined ranges for L and albedo
L_range = [1100, 1400]				# Range of solar constant
albedo_range = [0.15, 0.65]
albedo = albedo_range[0]

# Define plot type
# 'L': Plot of steady state temperature (K) vs. solar constant (W m-2)
# 'iterDown': Temperature (K) for each value of solar constant (W/m2) vs. iteration number
# 'iterUp': Temperature (K) for each value of solar constant (W/m2) vs. iteration number
plotType = 'L' # Options for type of plot: 'L', 'iterDown', 'iterUp'

# Initializing 
x, y = [],[]
L = L_range[1]
n = 30

while L >= L_range[0]:
	for i in range(0, n):
		temp = ((L * ( 1 - albedo ) ) / (4 * sigma)) ** 0.25
		albedo = albedo_gradient * temp + albedo_intercept
		albedo = min(albedo, albedo_range[1])
		albedo = max(albedo, albedo_range[0])
		if plotType is "iterDown":
			x.append(i)
			y.append(temp)
	if plotType == 'L':	
		x.append(L)
		y.append(temp)
	if plotType == 'iterDown':
		x.append(np.nan)
		y.append(np.nan)
	L -= 10

while L <= L_range[1]:
	for i in range(0, n):
		temp = ((L * ( 1 - albedo ) ) / (4 * sigma)) ** 0.25
		albedo = albedo_gradient * temp + albedo_intercept
		albedo = min(albedo, 0.65)
		albedo = max(albedo, 0.15)
		if plotType is "iterUp":
			x.append(i)
			y.append(temp)
	if plotType == 'L':	
		x.append(L)
		y.append(temp)
	if plotType == 'iterUp':
		x.append(np.nan)
		y.append(np.nan)
	L += 10


# Plot Graph
plt.plot( x, y )
plt.title('Albedo Hysterisis Loop')

if plotType == 'L':
	plt.xlabel('Solar Constant (W m-2)')
	plt.ylabel('Steady State Temperature (K)')
	plt.title('Albedo Hysterisis Loop')

if plotType == 'iterDown':
	plt.xlabel('Iteration Number')
	plt.ylabel('Temperature (K)')
    
if plotType == 'iterUp':
	plt.xlabel('Iteration Number')
	plt.ylabel('Temperature (K)')

plt.show( )
