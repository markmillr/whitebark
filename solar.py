import numpy as np
import time
import matplotlib.pyplot as plt

t0 = time.clock()

pi = np.pi

S = 1361 # Annual mean solar constant (W/m2)
albedo = 0.12 # Mean Bond albedo for Moon
Sabs = S * (1.0 - albedo)
P_moon = 29.53059*24.*3600. #Mean length of Moon SYNODIC day [s]

modelruntime = 2*P_moon
dt = 1.0 * 3600. # timestep in seconds (1.0*3600*24 = 1 earth day)
Nt = int(round(modelruntime/dt))   # Calculate number of timesteps (rounded and converted to integer)


t = np.zeros(Nt) # start the clock at zero









print(t, dt)

modeldays = 29.53059 # Earth days
modelruntime = modeldays*24.*3600. # in Earth seconds

Nt = int(round(modelruntime / dt)) # Number of timesteps


TWOPI = 2.0*np.pi

def hourAngle(t):

	return (TWOPI * t/day*24) % TWOPI

def cosSolarZenith(lat, dec, h):

	# Cosine of solar zenith angle
    x = np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(h)

    y = 0.5*(x + np.abs(x))

    return y

h_array = np.zeros(Nt)
y_array = np.zeros(Nt)

for n in range(0, Nt):
	t += dt
	h = hourAngle(t)
	c = cosSolarZenith(lat, dec, h)
	i = np.arccos(c)

	h_array[n] = t/day/24
	y_array[n] = h

	print(t/day*24, i, np.degrees(i))


# Stop timer for calculating execution time

t1 = time.clock()

print('\nExecution time:', t1-t0, 's')

plt.scatter(h_array, y_array)
plt.plot(h_array, y_array)
plt.xlabel('time in lunar hours')
plt.show()
