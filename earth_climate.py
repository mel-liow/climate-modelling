
##
# The future of Earth's climate
# A simulation of the near-term future of Earth's temperature as a function of 
# future CO2 emissions, climate sensitivity and the current cooling impact of
# industrial aerosols and short-lived greenhouse gases.
# There are two plots: 
# 1) Business as usual
# 2) World without humans 
##

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
import math

# Predefined constants
time_range = [1900, 2100]			# years
time_step = 1			# years
n = int((time_range[1] - time_range[0]) / time_step)			# number of data values

eq_CO2 = 280			# equilibrium CO2
i_CO2 = 290			# initial CO2

growth_rate = 0.0225			# rate that gives an atmospheric rise of 2.5ppm when the pCO2 = 400ppm
drawdown_rate = 0.01

watts_m2_2x = 4.0			# radiative forcing in Watts/m2 for doubling CO2
climate_sensitivity_2x = 5.9			# temperature change per Watts/m2 of forcing 
climate_sensitivity_watts_m2 = climate_sensitivity_2x / watts_m2_2x
t_response_time = 20

rf_aero_now = -1.5			# aerosol radiative forcing at present

# Initialising lists
p_CO2 = np.zeros(n)
p_CO2[0] = i_CO2
p_change = np.zeros(n)

rf_CO2 = np.zeros(n)
rf_masked = np.zeros(n)
rf_total = np.zeros(n)

T_eq = np.zeros(n)
T_transient = np.zeros(n)

years = np.arange(time_range[0], time_range[1])

# Loop through time
# A business as usual world
for i in range(1, len(years)):
	p_CO2[i] = eq_CO2 + (p_CO2[i-1] - eq_CO2) * (1 + growth_rate * time_step)
	p_change[i] = (p_CO2[i] - p_CO2[i-1]) / time_step
	rf_CO2[i] =  watts_m2_2x * np.log(p_CO2[i] / eq_CO2) / np.log(2)

iYear = years.tolist().index(2015)
aerosol_masking_factor = rf_aero_now / (p_CO2[iYear] - p_CO2[iYear - 1]) / time_step

for i in range(1, len(years)):
	rf_masked[i] = max(p_change[i] * aerosol_masking_factor, rf_aero_now)
	rf_total[i] = rf_masked[i] + rf_CO2[i]
	T_eq[i] = rf_total[i] * climate_sensitivity_watts_m2
	T_transient[i] = T_transient[i-1] + (T_eq[i] - T_transient[i-1]) * time_step/t_response_time

# A world without humans
p_CO2_rampdown = np.zeros(n)
p_CO2_rampdown[0:iYear] = p_CO2[0:iYear]

rf_CO2_rampdown = np.zeros(n)
rf_CO2_rampdown[0:iYear] = rf_CO2[0:iYear]
rf_total_rampdown = np.zeros(n)
rf_total_rampdown[0:iYear] = rf_total[0:iYear]

T_eq_rampdown = np.zeros(n)
T_eq_rampdown[0:iYear] = T_eq[0:iYear]
T_transient_rampdown = np.zeros(n)
T_transient_rampdown[0:iYear] = T_transient[0:iYear]

for i in range(iYear, len(years)):
	p_CO2_rampdown[i] = p_CO2_rampdown[i-1] + (eq_CO2 * 1.2 - p_CO2_rampdown[i-1]) * (drawdown_rate * time_step)
	rf_CO2_rampdown[i] = watts_m2_2x * np.log(p_CO2_rampdown[i]/eq_CO2) / np.log(2)
	rf_total_rampdown[i] = rf_CO2_rampdown[i]
	T_eq_rampdown[i] = rf_total_rampdown[i] * climate_sensitivity_watts_m2
	T_transient_rampdown[i] = T_transient_rampdown[i-1] + (T_eq_rampdown[i] - T_transient_rampdown[i-1]) * time_step / t_response_time


print(T_transient[iYear], T_transient_rampdown[iYear])

# Plot Graph
line1,=plt.plot(years,T_transient,label="Busines-as-usual")
line2,=plt.plot(years,T_transient_rampdown,'-',label="World without humans after 2015")
plt.legend(handler_map={line1:HandlerLine2D(numpoints=4)})

plt.title("Future Climate Model")
plt.xlabel("Year")
plt.ylabel("Transient Temperature (C)")
plt.grid()
plt.show()