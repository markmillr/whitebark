import numpy as np
import matplotlib.pyplot as plt

z_ = np.linspace(0, 1, 100)


rho_d = 1800# Deep layer density [kg m-3] Carrier et al 1991
rho_s = 1100# Deep layer density [kg m-3] Hayne et al 2013
rho_i = np.zeros(len(z_)) # density at any layer 
H = 0.06 # scale parameter from Hayne et al 2017

rho_i[:] = rho_d - (rho_d - rho_s)*np.exp(-z_/H)


fig, ax = plt.subplots()
ax.plot(rho_i, z_, color="#2242c7")
ax.set_ylim(np.amax(z_), np.amin(z_))
ax.set_xlabel("Soil density (kg/m3)")
ax.set_ylabel("Soil depth (m)")
ax.set_title("Lunar regolith thermal model soil density, scale parameter H = 0.06 m")
plt.show()
