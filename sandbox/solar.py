# solar.py

# Dependencies

import matplotlib.pyplot as plt    # for plotting
import numpy as np    # numpy arrays for calculations and matplotlib plots
from scipy import constants

# Intra-package references
import planets
import plotter

pi = planets.pi    # so we don't have to type np.pi all the time
sigma = planets.sigma

# Moon climate parameters

S_E = planets.Moon.S
albedo = planets.Moon.albedo
P_moon = planets.Moon.day
epsilon = planets.Moon.emissivity
q_c = planets.q_c
q_g = planets.Moon.Qb

Sabs = S_E * (1.0 - albedo)

# estimate mean equatorial soil heat flux from ballpark nighttime temp of 100 K
q_soil = epsilon*sigma*(100**4) - q_c - q_g
print('q_soil: {:.2f} W m-2'.format(q_soil))

# Model parameters
modelruntime = 2*P_moon
dt = 1.0 * 3600. # timestep in seconds (1.0*3600*24 = 1 earth day)
Nt = int(round(modelruntime/dt))   # Calculate number of timesteps (rounded and converted to integer)

# Create numpy arrays for calculations and plotting
t = np.zeros(Nt) # start the clock at zero
h = np.zeros(Nt) # hour angle
psi = np.zeros(Nt) # clipping function
solar = np.zeros(Nt) # array of insolation values to be calculated 
T_s = np.zeros(Nt) # time array of surface temperatures

# Create numpy array for lunar hour

lunarhour = np.zeros(Nt)

#print(int(P_moon/2*dt), Nt)
# Time loop, goes until modelruntime
for n in range(0, Nt):
	t[n] = n*dt    # Calculate the model time based on the current timestep number
	h[n] = 2*pi*(t[n]/P_moon) # Calculate hour angle
	psi[n] = 0.5 * ( np.cos(h[n]) + np.abs(np.cos(h[n]))) # Calculate clipping function
	solar[n] = Sabs * psi[n]  # Calculate solar flux

	#T_s[n] = ((1-albedo)*solar[n]/epsilon/sigma)**(1/4) # Calculate surface temperature given solar flux
	#T_s[n] = (((1-albedo)*solar[n] + q_c)/epsilon/sigma)**(1/4) # Calculate surface temperature given solar flux and cosmic flux
	#T_s[n] = (((1-albedo)*solar[n] + q_c + q_g)/epsilon/sigma)**(1/4) # Calculate surface temperature given solar, cosmic, and geothermal flux
	T_s[n] = (((1-albedo)*solar[n] + q_c + q_g + q_soil)/epsilon/sigma)**(1/4) # Calculate surface temperature given solar, cosmic, and geothermal flux

	lunarhour[n] = t[n]/P_moon*24 - 12  # Calculate lunar hour (for plotting only)




# Plot arrays using local plotter module

plot_title = 'Model lunar equatorial sfc temp with: S, q_c, q_g, q_soil'
plotter.plot(x=lunarhour, y=T_s, title = plot_title)

