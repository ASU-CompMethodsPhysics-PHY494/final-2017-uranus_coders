# ASU PHY 494 Project 2: Data for the TRAPPIST-1 system
# Copyright (c) Oliver Beckstein and Ian Kenney 2017
# Released under the MIT License

# values used from:
# https://solarsystem.nasa.gov/planets/
# http://www.astronomynotes.com/tables/tablesb.htm
# http://www.windows2universe.org/our_solar_system/planets_orbits_table.html

import numpy as np

yr2day = 365.25

#theta0 and distance information from https://omniweb.gsfc.nasa.gov/coho/helios/planet.html
#date of data is 01/01/2000

solar_system = {
    'mercury': {
        'mass':3.301e23,    # in kg
        'radius':2.4397e3,   # in km
        'eccentricity':0.20563593,
        'period':87.969,   # days
        'semi-major':0.3871 * 1e3,   # 1e-3 au
        'distance':0.466 * 1e3, # au
        'theta0':177.04, # degrees
        'phi0':-3.27 #degrees
    },
    'venus': {
        'mass':4.867e24,    # in kg
        'radius':6.0518e3,   # in km
        'eccentricity':0.00677672,
        'period':224.701,   # days
        'semi-major':0.7233 * 1e3,   # 1e-3 au
        'distance':0.720 * 1e3, # au
        'theta0':106.02, # degrees
        'phi0':3.28 #degrees
    },
    'earth':  {
        'mass':5.972e24,   # in kg
        'radius': 6.371e3,  # in km
        'eccentricity':0.01671123,
        'period':365.25,   # days
        'semi-major':1.000 * 1e3,   # 1e-3 au
        'distance':0.983 * 1e3, # au
        'theta0':23.92, # degrees
        'phi0':-3.28 #degrees

    },
    'mars': {
        'mass':6.4169e23,    # in kg
        'radius':3.3895e3,   # in km
        'eccentricity':0.0933941,
        'period':686.98,   # days
        'semi-major':1.5273 * 1e3,   # 1e-3 au
        'distance':1.391 * 1e3, # au
        'theta0':283.50, # degrees
        'phi0':-3.28 #degrees
    },
    'jupiter': {
        'mass':1.8981e27,    # in kg
        'radius':6.9911e4,   # in km
        'eccentricity':0.04838624,
        'period':11.862 * yr2day,   # days
        'semi-major':5.2028 * 1e3,   # 1e-3 au
        'distance':4.967 * 1e3, # au
        'theta0':320.54, # degrees
        'phi0':-3.28 #degrees
    },
    'saturn': {
        'mass':5.6832e26,    # in kg
        'radius':5.8232e4,   # in km
        'eccentricity':0.05386179,
        'period':29.457 * yr2day,   # days
        'semi-major':9.5388 * 1e3,   # 1e-3 au
        'distance':9.139 * 1e3, # au
        'theta0':330.11, # degrees
        'phi0':-3.28 #degrees
    },
    'uranus': {
        'mass':8.6810e25,    # in kg
        'radius':2.5362e4,   # in km
        'eccentricity':0.04725744,
        'period':84.011 * yr2day,   # days
        'semi-major':19.1914 * 1e3,   # 1e-3 au
        'distance':19.907 * 1e3, # au
        'theta0':240.65, # degrees
        'phi0':-3.28 #degrees
    },
    'neptune': {
        'mass':1.0241e26,    # in kg
        'radius':2.4622e4,   # in km
        'eccentricity':0.00859048,
        'period':164.79 * yr2day,   # days
        'semi-major':30.0611 * 1e3,   # 1e-3 au
        'distance':30.087 * 1e3, # au
        'theta0':227.96, # degrees
        'phi0':-3.28 #degrees
    },
    'sun': {
        'mass':1.989e30,     # in kg
        'radius': 695.8e3    # in km
    }
}

astronomical_unit = 149597870.700  # in km
au = astronomical_unit

# units
# - mass: in star masses (TRAPPIST-1)
# - period: in earth days (1 yr = 365.25 d)
# - radius: in earth radii
# - semi-major: in 1e-3 au (e.g., 11.11 -> 11.11e-3 au)
# - eccentricity: unit-less
planets = {
    'name': np.array(['Mercury', 
                      'Venus', 
                      'Earth', 
                      'Mars', 
                      'Jupiter', 
                      'Saturn', 
                      'Uranus', 
                      'Neptune']),
                      
    'mass': np.array([solar_system['mercury']['mass'] / solar_system['sun']['mass'],
                      solar_system['venus']['mass'] / solar_system['sun']['mass'],
                      solar_system['earth']['mass'] / solar_system['sun']['mass'],
                      solar_system['mars']['mass'] / solar_system['sun']['mass'],
                      solar_system['jupiter']['mass'] / solar_system['sun']['mass'],
                      solar_system['saturn']['mass'] / solar_system['sun']['mass'],
                      solar_system['uranus']['mass'] / solar_system['sun']['mass'],
                      solar_system['neptune']['mass'] / solar_system['sun']['mass']]),

    'period': np.array([solar_system['mercury']['period'],
                        solar_system['venus']['period'],
                        solar_system['earth']['period'],
                        solar_system['mars']['period'],
                        solar_system['jupiter']['period'],
                        solar_system['saturn']['period'],
                        solar_system['uranus']['period'],
                        solar_system['neptune']['period']]),   # Value at index 9 is supposed to be the velocity of the celestial object
                        
    'radius': np.array([solar_system['mercury']['radius'] / solar_system['earth']['radius'],
                        solar_system['venus']['radius'] / solar_system['earth']['radius'],
                        solar_system['earth']['radius'] / solar_system['earth']['radius'],
                        solar_system['mars']['radius'] / solar_system['earth']['radius'],
                        solar_system['jupiter']['radius'] / solar_system['earth']['radius'],
                        solar_system['saturn']['radius'] / solar_system['earth']['radius'],
                        solar_system['uranus']['radius'] / solar_system['earth']['radius'],
                        solar_system['neptune']['radius'] / solar_system['earth']['radius']]),

    'semi-major': np.array([solar_system['mercury']['semi-major'],
                            solar_system['venus']['semi-major'],
                            solar_system['earth']['semi-major'],
                            solar_system['mars']['semi-major'],
                            solar_system['jupiter']['semi-major'],
                            solar_system['saturn']['semi-major'],
                            solar_system['uranus']['semi-major'],
                            solar_system['neptune']['semi-major']]),

    'eccentricity': np.array([solar_system['mercury']['eccentricity'],
                              solar_system['venus']['eccentricity'],
                              solar_system['earth']['eccentricity'],
                              solar_system['mars']['eccentricity'],
                              solar_system['jupiter']['eccentricity'],
                              solar_system['saturn']['eccentricity'],
                              solar_system['uranus']['eccentricity'],
                              solar_system['neptune']['eccentricity']]),

    'distance': np.array([solar_system['mercury']['distance'],
                            solar_system['venus']['distance'],
                            solar_system['earth']['distance'],
                            solar_system['mars']['distance'],
                            solar_system['jupiter']['distance'],
                            solar_system['saturn']['distance'],
                            solar_system['uranus']['distance'],
                            solar_system['neptune']['distance']]),

    'theta0': np.array([solar_system['mercury']['theta0'],
                            solar_system['venus']['theta0'],
                            solar_system['earth']['theta0'],
                            solar_system['mars']['theta0'],
                            solar_system['jupiter']['theta0'],
                            solar_system['saturn']['theta0'],
                            solar_system['uranus']['theta0'],
                            solar_system['neptune']['theta0']]),

    'phi0': np.array([solar_system['mercury']['phi0'],
                            solar_system['venus']['phi0'],
                            solar_system['earth']['phi0'],
                            solar_system['mars']['phi0'],
                            solar_system['jupiter']['phi0'],
                            solar_system['saturn']['phi0'],
                            solar_system['uranus']['phi0'],
                            solar_system['neptune']['phi0']]),

    'rxy': np.array([solar_system['mercury']['distance'] * np.cos(solar_system['mercury']['phi0']),
                            solar_system['venus']['distance'] * np.cos(solar_system['venus']['phi0']),
                            solar_system['earth']['distance'] * np.cos(solar_system['earth']['phi0']),
                            solar_system['mars']['distance'] * np.cos(solar_system['mars']['phi0']),
                            solar_system['jupiter']['distance'] * np.cos(solar_system['jupiter']['phi0']),
                            solar_system['saturn']['distance'] * np.cos(solar_system['saturn']['phi0']),
                            solar_system['uranus']['distance'] * np.cos(solar_system['uranus']['phi0']),
                            solar_system['neptune']['distance'] * np.cos(solar_system['neptune']['phi0'])])
}

star = {
    'name': 'Sun',
    'mass': 1.0,            # in star masses
    'mass_solar': 1.0,   # in solar masses
    'radius_solar': 1.0,  # in solar radii
}


# mass of TRAPPIST-1 in solar masses
M_star_in_solar_mass = star['mass_solar']
# radius of TRAPPIST-1 in 1e-3 au
# 0.544182879201903  10^-3 au
star_radius_localunits = star['radius_solar'] * solar_system['sun']['radius'] / \
                         (astronomical_unit * 1e-3)


# in AU (solar system) units:
# GM = 4pi^2 AU^3/year^2
#G = 4*np.pi**2

# Gravitational constant in TRAPPIST-1 local units
# length = 1e-3 AU, time = d, M = M_star
G_local = 4*np.pi**2 * (1e3)**3 / (365.25)**2 * M_star_in_solar_mass
