##
# Temperature time-series of how the planetary temperature of a naked planet
# evloves over time to reach equilibrium where
# the incoming energy flux from the Sun equals outgoing energy flux from the planet
#
# Differential equation of heat content is given by
# dHeatContent/dt = Flux in - Flux out = L(1 - albedo) - εσT^4
# where HeatContent is a function of Temperature:
# HeatContent [J/m2] = Temperature [K] * HeatCapacity [J/m2 K]
##

import numpy as np
import matplotlib.pyplot as plt
import itertools as it
from collections import namedtuple


def temperature_timeseries(temp=0, heat=0, time=0, step=50, water_depth=4000, solar_constant=1350, albedo=0.3, epsilon=1, sigma=5.67E-8, specific_heat=4.2e3, density=997):

    heat_capacity = water_depth * density * specific_heat   # J/K
    flux_in = solar_constant * (1 - albedo) / 4             # W/m2
    dt = step * 365.25 * 3600 * 24                          # s

    while True:
        flux_out = epsilon * sigma * temp**4
        net_flux = flux_in - flux_out

        yield (time, temp, flux_out)

        heat += net_flux * dt
        time += step
        temp = heat / heat_capacity


# Plot temperature evolution
def plot_planet_evolution(data):
    time_list = [i for i, _, _ in data]
    temp_list = [t for _, t, _ in data]
    flux_out = [f for _, _, f in data]

    plt.plot(time_list, temp_list, label='Temperature')
    plt.title('Time-dependent Naked Planet Simulation')
    plt.xlabel('Time (years)')
    plt.ylabel('Temperature (K)')

    ax = plt.gca().twinx()
    ax.plot(time_list, flux_out, '-r', label='Outgoing heat flow')
    ax.set_ylabel('Flow (W/m^2)')

    plt.legend()


# Read input or return default value.
def read_input(user_input, default=None):
    data = input(user_input)
    if data:
        return int(data)
    elif default is not None:
        return default
    else:
        return read_input(user_input)


if __name__ == '__main__':
    steps = read_input('# of half centuries: ', 25)

    gen = temperature_timeseries()
    data = []
    for i in range(0, steps):
        data.append(list(next(gen)))

    plot_planet_evolution(data)
    plt.show()
