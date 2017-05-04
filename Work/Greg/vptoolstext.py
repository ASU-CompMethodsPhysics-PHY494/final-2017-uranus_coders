from vpython import *
from PIL import ImageGrab
from subprocess import call
import numpy as np
import math
import main
import parameters
import matplotlib
import matplotlib.pyplot as plt

t_max = int(2*parameters.solar_system['jupiter']['period'])
ct = int(1*parameters.solar_system['jupiter']['period'])
dt = 1
scene2 = canvas(title="Solar System")
scene3 = canvas(title="Timer")

def simulate(t_max=t_max, dt=dt, ct=ct):
	nsteps = int(t_max/dt)
	time, r, v, coPos = main.integrate_orbits(dt=dt, t_max=t_max, collision_time=ct)
	#rA[i,j,k]
	#i = timestep
	#j = planet
	#k = x,y,z
	rA = np.insert(r[:], 2, 0, axis=2)
	rA = np.insert(rA[:], 8, 0, axis=1)
	vA = np.insert(v[:], 2, 0, axis=2)
	rCO = np.insert(coPos[:], 2, 0, axis=1)
	rCOs = np.zeros_like(rA[:,0,:])

	for i in range(nsteps):
		if i < ct:
			rA[i,8] = rCO[i]
		else:
			rA[i,8] = rA[i,4,:]


	sun = sphere(pos=vector(0,0,0), color=color.orange, make_trail=True, canvas=scene2, radius=50)
	mercury = sphere(pos=vector(rA[0,0,0],rA[0,0,1],rA[0,0,2]), color=color.gray(0.5), make_trail=True, canvas=scene2, radius=50)
	venus = sphere(pos=vector(rA[0,1,0],rA[0,1,1],rA[0,1,2]), color=color.yellow, make_trail=True, canvas=scene2, radius=50)
	earth = sphere(pos=vector(rA[0,2,0],rA[0,2,1],rA[0,2,2]), color=color.blue, make_trail=True, canvas=scene2, radius=50)
	mars = sphere(pos=vector(rA[0,3,0],rA[0,3,1],rA[0,3,2]), color=color.red, make_trail=True, canvas=scene2, radius=50)
	jupiter = sphere(pos=vector(rA[0,4,0],rA[0,4,1],rA[0,4,2]), color=color.orange, make_trail=True, canvas=scene2, radius=50)
	saturn = sphere(pos=vector(rA[0,5,0],rA[0,5,1],rA[0,5,2]), color=color.green, make_trail=True, canvas=scene2, radius=50)
	uranus = sphere(pos=vector(rA[0,6,0],rA[0,6,1],rA[0,6,2]), color=color.cyan, make_trail=True, canvas=scene2, radius=50)
	neptune = sphere(pos=vector(rA[0,7,0],rA[0,7,1],rA[0,7,2]), color=color.white, make_trail=True, canvas=scene2, radius=50)
	CO = sphere(pos=vector(rA[0,8,0],rA[0,8,1],rA[0,8,2]), color=color.magenta, make_trail=True, canvas=scene2, radius=50)
	planets = (mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, CO)

	return rA, vA, planets, nsteps

"""def timer(nsteps):
	k = 0
	day = np.zeros(nsteps)
	year = np.zeros_like(day)
	jyear = np.zeros_like(day)
	cyear = np.zeros_like(day)
	while k < nsteps:
		day[k] = str(k)
		year[k] = str(int(100*k/365)/100)
		jyear[k] = str(int(100*k/(365*11.86))/100)
		cyear[k] = str(int((2*(100*parameters.solar_system['jupiter']['period']-k))/365)/100)
		k += 1
	return day, year, jyear, cyear"""

def solplot(t_max=t_max, dt=dt, ct=ct, r8=10):
	rA, vA, planets, nsteps = simulate(t_max=t_max, dt=dt, ct=ct)
	i = 0
	fnum = 0
	"""day, year, jyear, cyear = timer(nsteps)
	text(pos=vector(-15, 15,0), text='CO mass: ' + str(parameters.SCALE_FACTOR) + "Jupiter Masses", canvas=scene3)
	text(pos=vector(-15, 5,0), text='CO velocity: ' + str(parameters.celestial_object['velocity_magnitude']) + ' AU/day', canvas=scene3)
	text(pos=vector(-15, -5,0), text='Collision time: ' + str(ct) + ' Jovian Years', canvas=scene3)
	text(pos=vector(-15, -15,0), text='Simulation time: ' + str(t_max) + ' Jovian Years',canvas=scene3)"""
	while i < nsteps:
		rate(r8)
		j = 0
		while j < len(planets):
			planets[j].pos = vector(rA[i,j,0], rA[i,j,1], rA[i,j,2])
			j += 1
		"""text(pos=vector(15, 15,0), text='Day: ' + str(day[i]), canvas=scene3)
		text(pos=vector(15, 5,0), text='Year: ' + str(year[i]), canvas=scene3)
		text(pos=vector(15, -5,0), text='Jovian Year: ' + str(jyear[i]), canvas=scene3)
		text(pos=vector(15, -15,0), text='Years to impact: ' + str(cyear[i]), canvas=scene3)"""
		i += 1	