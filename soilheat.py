
#   1-D thermal conduction model of lunar regolith. Outputs lunar surface temperature.
#
#   Author: Mark A. Miller (github.com/markmillr)
#
#   License: MIT (./LICENSE.txt)
#
#   Adapted from numerical thermal model documented in: 
#
#   Hayne, P. O., Bandfield, J. L., 
#   Siegler, M. A., Vasavada, A. R., Ghent, R. R., Williams, J.-P., ... Paige, D. A. (2017). 
#   Global regolith thermophysical properties of the Moon from the Diviner Lunar Radiometer 
#   Experiment. Journal of Geophysical Research: Planets, 122, 2371–2400. 
#   https://doi.org/10.1002/ 2017JE005387
#
#   Current features:
#
#   Uniform grid spacing dz
#   Nonlinear (temperature-driven) soil thermodynamic properties (conductivity, thermal diffusivity)
#
#   TODO:
#
#   Implement geometrically increasing grid spacing dz (requires different numerical solution)
#   Access planetary parameters as an import from the planets.py database
#   Create Planet and Model classes for reusability
#   
#   Long-range: Create GUI


import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import scipy.optimize
import sympy
import math

# Set constants
pi = np.pi
sigma = spc.sigma

# Moon parameters

S = 1361. # Annual mean solar constant (W/m2)
albedo = 0.12 # Mean Bond albedo for Moon
epsilon = 0.98 # T7 emissivity (Vasavada et al.) 
Sabs = S * (1.0 - albedo)
P_moon = 29.53059*24.*3600. #Mean length of Moon SYNODIC day [s]
a = 1e-08   # Value for lunar regolith thermal diffusivity provided by Hayne et al 2017
#a = 1e-5 # try unrealistically high thermal diffusivity
k = 0.6e-03     # Equatorial thermal conductivity in Wm-2K-1 (summarized by Hayne et al 2017 Table A2)
q_g = 11e-3 # Mean lunar geothermal surface flux [W m^-2] (Siegler and Smrekar (2014)
q_c = 3.13e-6 # Cosmic background radiation flux [W m^-2] (Fixsen, 2009)

# Nonlinear Moon parameters (from Hayne et al 2017)

rho_d = 1800 # Deep layer density [kg m-3] Carrier et al 1991
rho_s = 1100 # Surface layer density [kg m-3] Hayne et al 2013
K_d = 3.4e-3 # Deep layer thermal conductivity [W m-1 K-1] Hayne et al 2017
K_s = 7.4e-4 # Surface layer thermal conductivity [W m-1 K-1] Hayne et al 2017
c0 = -3.6125 # Coefficients for specific heat capacity function [J kg-1 K-1] Hayne et al 2017
c1 = 2.7431
c2 = 2.3616e-3
c3 = -1.2340e-5
c4 = 8.9093e-9
chi = 2.7 # Radiative conductivity parameter from Hayne et al 2017

# Model parameters

T_sfc = 250.0 # Initial temperature [K] at top boundary above the zeroth layer
T_bot = 250.0 # Initial temperature [K] at bottom boundary below the last layer

modelruntime = 2*24*P_moon/24
dt = 3600#1.0*3600 # timestep in seconds (1.0*3600*24 = 1 earth day)
Nt = int(round(modelruntime/dt))   # Calculate number of timesteps (rounded and converted to integer)

D = 1.0     # Thickness (depth) of entire slab/system [m] 
dz = 0.3    # Thickness of each layer (dist. btwn gridpoints) [m]
Nz = int(round(D/dz))   # Calculate number of layers (rounded and converted to integer)

F = a*dt/(dz**2)    # Mesh Fourier number - term in heat conduction equation

Ts = np.float() # surface temperature variable used in top boundary condition 

# Create np arrays for solar calculations and plotting
t = np.zeros(Nt) # start the clock at zero
lt = np.zeros(Nt) # local time in local hours
h = np.zeros(Nt) # hour angle
psi = np.zeros(Nt) # clipping function
solar = np.zeros(Nt) # time array of insolation values to be calculated 
surfaceflux = np.zeros(Nt) # time array of solar + cosmic 
Ts_array = np.zeros(Nt) # time array of surface temperatures

# Create np arrays for soil thermodynamic properties
rho = np.zeros(Nz+1) # array of densities
K_c = np.zeros(Nz+1) # array of thermal conductivities
K = np.zeros(Nz+1) # array 
c = np.zeros(Nz+1) # array of specific heat capacities

# Create np arrays for current soil temperature (T_n) and new soil temperature at timestep n+1 (T)

T_n = np.zeros(Nz+1)
T = np.zeros(Nz+1)

# Create np array for soil depth (for plotting)

z = np.zeros(Nz+1)

# Set initial condition T(z,0) and fill in values for depth array (for plotting)

for i in range(0, Nz+1):
    if i == 0:
        T_n[i] = T_sfc # Ballpark equatorial max temp [K]
    elif i == Nz:
        T_n[i] = T_bot # Ballpark day/night mean temp for bottom layer [K]
    else:
        T_n[i] = T_bot # Unrealistic: set all remaining layers to T_bottom
    z[i] = i*dz

# Time-stepping loop

for n in range(0, Nt):
    # Calculate surface flux from above
    t[n] = n*dt # Calculate model time elapsed
    lt[n] = t[n]/P_moon*24 - 12 # Calculate local time in local hours
    h[n] = 2*pi*(t[n]/P_moon) # Calculate hour angle
    psi[n] = 0.5 * ( np.cos(h[n]) + np.abs(np.cos(h[n]))) # Calculate clipping function
    solar[n] = Sabs * psi[n]  # Calculate solar flux
    surfaceflux[n] = solar[n] + q_c # Calculate total flux of energy downward onto surface - Probably useless

    # Calculate soil thermal properties at this timestep
    
    # Scalar version
    
    for i in range(0, Nz+1):
        rho[i] = rho_d - (rho_d - rho_s)*math.exp(-i/D)
        K_c[i] = K_d - (K_d-K_s)*(rho_d - rho[i])/(rho_d - rho_s)
        K[i] = K_c[i] * (1 + chi* ( (T_n[i]/350)**3 )   )
        c[i] = c0 + c1*T_n[i] + c2*(T_n[i]**2) + c3*(T_n[i]**3) + c4*(T_n[i]**4)
    
    # Set top boundary condition (Robin)

    Ts = T_n[0]

    T[0] = Ts

    #Coefficients of top boundary condition in polynomial form
    Kc0 = K_c[0]
    a0 = epsilon*sigma + 3/2/dz*Kc0*chi/(350**3)
    a1 = Kc0*chi/2/dz/(350**3)*(4*T_n[1] - T_n[2])
    a2 = 0
    a3 = -3*Kc0/2/dz
    a4 = -3*Kc0/2/dz*(4*T_n[1]-T_n[2]) - surfaceflux[n]

    def f(Ts):
            return a0*Ts**4 + a1*Ts**3 + a2*Ts**2 + a3*Ts + a4
        
    try:
        Ts = abs(scipy.optimize.newton(f, Ts))
        method = "Newton"

    except:
        print("root solver not happy")

    # Set bottom boundary condition
    T[Nz] = T_bot
    

    # Set top layer temp
    T[0] = Ts

    # Compute u at inner mesh points (vectorized, nonlinear diffusivity)
    
    T[1:Nz-1] = T_n[1:Nz-1]+ dt/rho[1:Nz-1]/c[1:Nz-1] * (K[1:Nz-1]/dz*(T[2:Nz] - T[1:Nz-1]) - K[0:Nz-2]/dz*(T[1:Nz-1] - T[0:Nz-2]))


    # Switch variables before next timestep
    #T_n, T = T, T_n
    T_n = T

    # Populate plotting array

    Ts_array[n] = Ts
    print(n*dt/P_moon*24, surfaceflux[n], Ts, T_n[0], T_n[1], T_n[2], T_n[3], method)

    
fig, ax = plt.subplots()
ax.scatter(lt, Ts_array)
ax.plot(lt, Ts_array)
ax.set_xlim(0,24)
ax.set_ylim(-20,400)
ax.set_xticks(np.arange(0, 25, step=6))
ax.set_xlabel("Time (lunar hours)")
ax.set_ylabel("Surface temperature (K)")
ax.annotate("Max temp: {:.0f} K".format(np.amax(Ts_array)), xy=(17,350))
ax.annotate("Min temp: {:.0f} K".format(np.amin(Ts_array)), xy=(17,325))
ax.annotate("dz: {} m".format(dz), xy=(1, 350))
ax.set_title("Model lunar equatorial sfc temp, dynamic soil thermal model included")
plt.show()






