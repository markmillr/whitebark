# Generates various plots showing dependence of model calculated conductivity k on temperature and density

# Refers to heat1d thermal model for lunar regolith published by 

#   Hayne, P. O., Bandfield, J. L., 
#   Siegler, M. A., Vasavada, A. R., Ghent, R. R., Williams, J.-P., ... Paige, D. A. (2017). 
#   Global regolith thermophysical properties of the Moon from the Diviner Lunar Radiometer 
#   Experiment. Journal of Geophysical Research: Planets, 122, 2371–2400. 
#   https://doi.org/10.1002/ 2017JE005387


import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap

z_ = np.linspace(0., 1., 100)
T_ = np.linspace(0, 400, 400)

T = 250

rho_d = 1800# Deep layer density [kg m-3] Carrier et al 1991
rho_s = 1100# Deep layer density [kg m-3] Hayne et al 2013
rho_i = np.zeros(len(z_))
chi = 2.7 # Radiative conductivity parameter Hayne et al 2017 and Vasavada et al 2012
k_d = 3.4e-3 # Deep layer thermal conductivity [W m-1 K-1] Hayne et al 2017
k_s = 7.4e-4 # Surface layer thermal conductivity [W m-1 K-1] Hayne et al 2017


H = 0.06 # scale parameter from Hayne et al 2017

# Density at each depth
rho_i[:] = rho_d - (rho_d - rho_s)*np.exp(-z_/H)

# Calculation of k as a function of T at constant surface density:

k_Ci = np.zeros(len(T_))
k_i = np.zeros(len(T_))

# Surface contact conductivity:
k_Ci[:] = k_d - (k_d - k_s)*( (rho_d - rho_s) / (rho_d - rho_s) )

# Radiation and contact conductivity:
k_i[:] = k_Ci[:]*(1 - chi*(T_[:]/350)**3)

fig, ax = plt.subplots()
ax.plot(T_, k_i, color="#2242c7")
ax.set_xlim(0, 400)
ax.axhline(y=0, linestyle="--", color="grey", linewidth="1.25")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Radiation and contact thermal conductivity (W/m/K)")
ax.set_title("\n".join(wrap("Lunar regolith thermal model dependence of thermal conductivity k on T at constant sfc density", 60)))
plt.show()

# Calculation of k as a function of rho at constant temperature 250 K:

k_Ci = np.zeros(len(rho_i))
k_i = np.zeros(len(rho_i))

# Surface contact conductivity:
k_Ci[:] = k_d - (k_d - k_s)*( (rho_d - rho_i[:]) / (rho_d - rho_s) )

# Radiation and contact conductivity:
k_i[:] = k_Ci[:]*(1 - chi*(T/350)**3)

fig, ax = plt.subplots()
ax.plot(rho_i, k_i, color="#2242c7")
ax.set_xlim(rho_s, rho_d)
ax.axhline(y=0, linestyle="--", color="grey", linewidth="1.25")
ax.set_xlabel("Soil density (kg/m3)")
ax.set_ylabel("Radiation and contact thermal conductivity (W/m/K)")
ax.set_title("\n".join(wrap("Lunar regolith thermal model dependence of thermal conductivity k on density at constant temperature 250 K", 60)))
plt.show()


"""
fig, ax = plt.subplots()
ax.plot(rho_i, z_, color="#2242c7")
ax.set_ylim(np.amax(z_), np.amin(z_))
ax.set_xlabel("Soil density (kg/m3)")
ax.set_ylabel("Soil depth (m)")
ax.set_title("Lunar regolith thermal model soil density, scale parameter H = 0.06 m")
plt.show()
"""
