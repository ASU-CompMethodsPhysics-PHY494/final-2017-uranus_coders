# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 13:35:42 2017

@author: Vishal
"""
import math
import numpy as np
import parameters

# ALL PLANETS WILL START FROM -1 * X-AXIS WITH INITIAL VELOCITY ONLY IN Y-AXIS [0, v_0]
# PLANETS TRAVELING CLOCKWISE

# \description Initial position r_max = r_max = semi_major_a * (1 + eccentricity). 
# \return [-r_max, 0]

#I commented out your initial condition functions and made my own

"""def initial_position(a, ecc):
    if(ecc == 0.0):
        # TODO - Determine appropriate initial position
        return np.array(0, parameters.planets['semi-major'][8])  # Object 9 is celestial object. Value is initial distance away from sun instead of semi-major axis value.
        
    r_max = a * (1 + ecc)
    return np.array([-r_max, 0.0])

# \description Initial velocity at aphelion
# \return [0, v_0]
def initial_velocity(a, ecc):
    if(ecc == 0.0):
        # TODO - Determine appropriate initial velocity
        return np.array(0, parameters.planets['period'][8])  # Object 9 is celestial object. Value is initial velocity instead of period
    v2 = parameters.G_local * parameters.star['mass'] * (1 - ecc) / (a * (1 + ecc))
    v_ap = math.sqrt(v2)
    return np.array([0, v_ap])
"""

def initial_position(angle, distance):
    x = np.deg2rad(angle)
    return distance * np.array([np.cos(x), np.sin(x)])

def initial_velocity(angle, distance, period):
    x = np.deg2rad(angle)
    v = 2*np.pi*distance/period
    return v * np.array([np.sin(x), -np.cos(x)])

# \description Calculate the force between 2 objects
# \return Force    
def F_gravity(r, m, M):
    """Force due to gravity between two masses.

    Parameters
    ----------
    r : array
      distance vector (x, y)
    m, M : float
      masses of the two bodies

    Returns
    -------
    array
       force due to gravity (along r)
    """
    rr = np.sum(r*r)
    rhat = r/np.sqrt(rr)
    return -parameters.G_local * m * M / rr * rhat
    
# \description Calculate net force on each planet based on current positions
# \return Array of net forces    
def F_planets(positions):
    length = len(parameters.planets['name'])
    retVal = np.zeros((length, 2))
    force = np.zeros((length, length, 2))    # 2-D array of individual forces
    
    # Fill in lower half of the array since upper half is equal value but in opposite direction
    for i in range(length):
        for j in range(i + 1):
            if(i == j): # Force on star
                force[i][j] = F_gravity(positions[i], parameters.planets['mass'][i], parameters.star['mass'])
            else:
                force[i][j] = F_gravity(positions[i] - positions[j], parameters.planets['mass'][i], parameters.planets['mass'][j])
                
    # Fill in the upper-half of array
    for i in range(length):
        for j in range(i + 1, length):
            force[i][j] = -1 * force[j][i]
            
    # Calculate net force
    for i in range(length):
        netForce = np.array([0.0, 0.0])
        for j in range(length):
            netForce += force[i][j]
        
        retVal[i] = netForce
        
    return retVal
        
# \description t_max in days
# \return position and velocity values of planets. [timeStep, numPlanets, positions/velocity]
def integrate_orbits(dt=1, t_max=1000):
    nsteps = int(t_max/dt)
    time = dt * np.arange(nsteps)

    # shape = (step, planet, x/y)
    r = np.zeros((nsteps, len(parameters.planets['name']), 2))
    v = np.zeros_like(r)
    
    length = len(parameters.planets['name'])
    # Set initial conditions
    for i in range(length):
        r[0, i, :] = initial_position(parameters.planets['semi-major'][i], parameters.planets['eccentricity'][i])
        v[0, i, :] = initial_velocity(parameters.planets['semi-major'][i], parameters.planets['eccentricity'][i], parameters.planets['period'][i])

    masses = np.zeros((length, 1))
    for i in range(length):
        masses[i] = [parameters.planets['mass'][i]]
    # start force evaluation for first step
    Ft = F_planets(r[0])
    for i in range(nsteps-1):
        vhalf = v[i] + 0.5*dt * Ft/masses
        r[i+1, :] = r[i] + dt * vhalf
        Ftdt = F_planets(r[i+1])
        v[i+1, :] = vhalf + 0.5*dt * Ftdt/masses
        # new force becomes old force
        Ft = Ftdt    
    
    return time, r, v

if __name__ == "__main__":    
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.style.use("ggplot")
    
    t, r, v = integrate_orbits(dt = .05, t_max = 1.5)
    rMercury = r[:,0]
    rVenus = r[:,1]
    rEarth = r[:,2]
    rMars = r[:,3]
    rJupiter = r[:,4]
    rSaturn = r[:,5]
    rUranus = r[:,6]
    rNeptune = r[:,7]
    #rCO = r[:,8]
    
    ax = plt.subplot()
    ax.plot(0, 0, "ro", label="Sun")
    ax.plot(rMercury[:,0], rMercury[:,1], label="Planet Mercury")
    ax.plot(rVenus[:,0], rVenus[:,1], label="Planet Venus")
    ax.plot(rEarth[:,0], rEarth[:,1], label="Planet Earth")
    ax.plot(rMars[:,0], rMars[:,1], label="Planet Mars")
    ax.plot(rJupiter[:,0], rJupiter[:,1], label="Planet Jupiter")
    ax.plot(rSaturn[:,0], rSaturn[:,1], label="Planet Saturn")
    ax.plot(rUranus[:,0], rUranus[:,1], label="Planet Uranus")
    ax.plot(rNeptune[:,0], rNeptune[:,1], label="Planet Neptune")
    #ax.plot(rCO[:,0], rCO[:,1], label="Celestial Object")
    ax.legend(loc="best")
    ax.set_xlabel(r"$x$ (AU)")
    ax.set_ylabel(r"$y$ (AU)")
    ax.set_title("Orbits for 1.5 days")
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
    ax = plt.subplot()
    t, r, v = integrate_orbits(dt = 1, t_max = 100000)
    star = np.zeros(len(r))
    rMercury = r[:,0]
    rVenus = r[:,1]
    rEarth = r[:,2]
    rMars = r[:,3]
    rJupiter = r[:,4]
    rSaturn = r[:,5]
    rUranus = r[:,6]
    rNeptune = r[:,7]
    #rCO = r[:,8]

    ax.plot(0, 0, "ro", label="Sun")
    ax.plot(rMercury[:,0], rMercury[:,1], label="Planet Mercury")
    ax.plot(rVenus[:,0], rVenus[:,1], label="Planet Venus")
    ax.plot(rEarth[:,0], rEarth[:,1], label="Planet Earth")
    ax.plot(rMars[:,0], rMars[:,1], label="Planet Mars")
    ax.plot(rJupiter[:,0], rJupiter[:,1], label="Planet Jupiter")
    ax.plot(rSaturn[:,0], rSaturn[:,1], label="Planet Saturn")
    ax.plot(rUranus[:,0], rUranus[:,1], label="Planet Uranus")
    ax.plot(rNeptune[:,0], rNeptune[:,1], label="Planet Neptune")
    #ax.plot(rCO[:,0], rCO[:,1], label="Celestial Object")
    ax.legend(loc="best")
    ax.set_xlabel(r"$x$ (AU)")
    ax.set_ylabel(r"$y$ (AU)")
    ax.set_title("Solar System Orbits")     
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    