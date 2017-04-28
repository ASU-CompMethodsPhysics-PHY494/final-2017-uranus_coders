from vpython import *
import numpy as np
import math
import main_mine as main
import matplotlib
import matplotlib.pyplot as plt

def simulate(t_max = 433, dt = 0.1):
	nsteps = int(t_max/dt)
	time, r, v = main.integrate_orbits(dt=dt, t_max=t_max)
	#rA[i,j,k]
	#i = timestep
	#j = planet
	#k = x,y,z
	rA = np.insert(r[:], 2, 0, axis=2)
	vA = np.insert(v[:], 2, 0, axis=2)

	scene2 = canvas(title="Solar System")

	sun = sphere(pos=vector(0,0,0), color=color.orange, make_trail=True, canvas=scene2, radius=50)
	mercury = sphere(pos=vector(rA[0,0,0],rA[0,0,1],rA[0,0,2]), color=color.gray(0.5), make_trail=True, canvas=scene2, radius=50)
	venus = sphere(pos=vector(rA[0,1,0],rA[0,1,1],rA[0,1,2]), color=color.yellow, make_trail=True, canvas=scene2, radius=50)
	earth = sphere(pos=vector(rA[0,2,0],rA[0,2,1],rA[0,2,2]), color=color.blue, make_trail=True, canvas=scene2, radius=50)
	mars = sphere(pos=vector(rA[0,3,0],rA[0,3,1],rA[0,3,2]), color=color.red, make_trail=True, canvas=scene2, radius=50)
	jupiter = sphere(pos=vector(rA[0,4,0],rA[0,4,1],rA[0,4,2]), color=color.orange, make_trail=True, canvas=scene2, radius=50)
	saturn = sphere(pos=vector(rA[0,5,0],rA[0,5,1],rA[0,5,2]), color=color.green, make_trail=True, canvas=scene2, radius=50)
	uranus = sphere(pos=vector(rA[0,6,0],rA[0,6,1],rA[0,6,2]), color=color.cyan, make_trail=True, canvas=scene2, radius=50)
	neptune = sphere(pos=vector(rA[0,7,0],rA[0,7,1],rA[0,7,2]), color=color.white, make_trail=True, canvas=scene2, radius=50)
	planets = (mercury, venus, earth, mars, jupiter, saturn, uranus, neptune)

	return rA, vA, planets, nsteps

def solplot(t_max = 433, dt = 0.1, r8 = 200):
	rA, vA, planets, nsteps = simulate(t_max = t_max, dt = dt)
	i = 0
	while i < nsteps-1:
		rate(r8)
		j = 0
		while j < len(planets):
		    planets[j].pos = vector(rA[i,j,0], rA[i,j,1], rA[i,j,2])
		    j += 1
		i += 1	