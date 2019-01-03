# solar.py

import matplotlib.pyplot as plt    # for plotting
import numpy as np    # numpy arrays for calculations and matplotlib plots
from scipy import constants


pi = np.pi    # so we don't have to type np.pi all the time
sigma = constants.sigma

# Moon climate parameters

S_E = 1361.1 # Solar constant at Earth TOA, +/- 0.5 [W m-2] (Gueymard, 2018)
albedo = 0.115 # Mean bond albedo value (ranges from 0.07 to 0.16) from Vasavada et al. (2012)
P_moon = 29.53059*24.*3600. #Mean length of Moon SYNODIC day [s]
epsilon = 0.98 	# T7 emissivity (Vasavada et al., 2012)
q_c = 3.13e-6 # Cosmic background radiation flux [W m^-2] (Fixsen, 2009)
q_g = 11e-3 # Mean lunar geothermal surface flux [W m^-2] (Siegler and Smrekar (2014)

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




# Plot solar and t arrays with matplotlib

"""
fig, ax = plt.subplots()
ax.scatter(lunarhour, T_s)
ax.plot(lunarhour, T_s)
ax.set_xlim(0,24)
ax.set_ylim(-20,400)
ax.set_xticks(np.arange(0, 25, step=6))
ax.set_xlabel('Time (lunar hours)')
ax.set_ylabel('Surface temperature (K)')
ax.annotate('Max temp: {:.0f} K'.format(np.amax(T_s)), xy=(17,350))
ax.annotate('Min temp: {:.0f} K'.format(np.amin(T_s)), xy=(17,325))
ax.set_title('Model lunar equatorial sfc temperature, solar heating only')
plt.show()
"""

fig, ax = plt.subplots()
ax.scatter(lunarhour, T_s)
ax.plot(lunarhour, T_s)
ax.set_xlim(0,24)
ax.set_ylim(-20,400)
ax.set_xticks(np.arange(0, 25, step=6))
ax.set_xlabel('Time (lunar hours)')
ax.set_ylabel('Surface temperature (K)')
ax.annotate('Max temp: {:.0f} K'.format(np.amax(T_s)), xy=(17,350))
ax.annotate('Min temp: {:.0f} K'.format(np.amin(T_s)), xy=(17,325))
#ax.annotate('Min temp: {:.1f} K'.format(np.amin(T_s)), xy=(17,325))
#ax.set_title('Model lunar equatorial sfc temperature, solar + cosmic heating')
#ax.set_title('Model lunar equatorial sfc temp, solar, cosmic, geothermal heating')
ax.set_title('Model lunar equatorial sfc temp with: S, q_c, q_g, q_soil')
plt.show()


"""
ax.plot(lunarhour, solar)
ax.xlabel('Time (lunar hours)')
ax.ylabel('Solar irradiance (W/m2)')
ax.xlim(0,24)
ax.xticks(np.arange(0, 25, step=6))
plt.show()


plt.scatter(lunarhour, solar)
plt.plot(lunarhour, solar)
plt.xlabel('Time (lunar hours)')
plt.ylabel('Solar irradiance (W/m2)')
plt.xlim(0,24)
plt.xticks(np.arange(0, 25, step=6))
plt.show()
"""

"""
plt.scatter(lunarhour, T_s)
plt.plot(lunarhour, T_s)
plt.xlabel('Time (lunar hours)')
plt.ylabel('Surface temperature (K)')
plt.xlim(0,24)
plt.ylim(-20, 400)
plt.xticks(np.arange(0, 25, step=6))
plt.show()
"""
