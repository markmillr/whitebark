# model.py

# Defines model parameters

import numpy as np
import planets

class Model(object):


	# Constructor method

	def __init__(self, planet=planets.Moon, ndays=1):

		# Initialize
		self.planet = planet
		self.Sabs = self.planet.S * (1.0 - self.planet.albedo)
		
		self.modelruntime = ndays*planet.day

		self.t = 0.
		self.dt = 1.0*3600. # timestep length in seconds
		self.Nt = int(round(self.modelruntime/self.dt)) # number of timesteps (rounded integer)

		self.h = 0. # hour angle
		self.psi = 0. # clipping function based on cosine of zenith angle
		self.q_solar = 0. # insolation value at point

		# Initialize arrays for output temperatures and local times

		#self.t = np.zeros(self.Nt) # Sidereal time
		self.lt = np.zeros(self.Nt) # local time [local hr]
		self.T = np.zeros(self.Nt) # output array of temps for plotting

	def run(self):

		for n in range(0, self.Nt):
			self.advance()
			self.lt[n] = self.t/self.planet.day*24.0 

	def advance(self):
		self.t += self.dt


	def solarSurfaceFlux(self):
		q_solar = self.Sabs*psi(self.h)

	def psi(self):
		self.h

modelrun = Model()
print(modelrun.modelruntime)


