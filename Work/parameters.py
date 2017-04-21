# ASU PHY 494 Project 2: Data for the TRAPPIST-1 system
# Copyright (c) Oliver Beckstein and Ian Kenney 2017
# Released under the MIT License

# values used from:
# https://solarsystem.nasa.gov/planets/
# http://www.astronomynotes.com/tables/tablesb.htm

import numpy as np

yr2day = 365.25

solar_system = {
    'mercury': {
        'mass':3.301e23,    # in kg
        'radius':2.4397e3,   # in km
        'eccentricity':0.20563593,
        'period':87.969,   # days
        'semi-major':   # au
    },
    'venus': {
        'mass':4.867e24,    # in kg
        'radius':6.0518e3,   # in km
        'eccentricity':0.00677672,
        'period':224.701,   # days
        'semi-major':   # au
    },
    'earth':  {
        'mass':5.972e24,   # in kg
        'radius': 6.371e3,  # in km
        'eccentricity':0.01671123,
        'period':365.25,   # days
        'semi-major':   # au
    },
    'mars': {
        'mass':6.4169e23,    # in kg
        'radius':3.3895e3,   # in km
        'eccentricity':0.0933941,
        'period':686.98,   # days
        'semi-major':   # au
    },
    'jupiter': {
        'mass':1.8981e27,    # in kg
        'radius':6.9911e4,   # in km
        'eccentricity':0.04838624,
        'period':11.862 * yr2day,   # days
        'semi-major':   # au
    },
    'saturn': {
        'mass':5.6832e26,    # in kg
        'radius':5.8232e4,   # in km
        'eccentricity':0.05386179,
        'period':29.457 * yr2day,   # days
        'semi-major':   # au
    },
    'uranus': {
        'mass':8.6810e25,    # in kg
        'radius':2.5362e4,   # in km
        'eccentricity':0.04725744,
        'period':84.011 * yr2day,   # days
        'semi-major':   # au
    },
    'neptune': {
        'mass':1.0241e26,    # in kg
        'radius':2.4622e4,   # in km
        'eccentricity':0.00859048,
        'period':164.79 * yr2day   # days
        'semi-major':   # au
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
    'name': np.array(['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 
                      'Uranus', 'Neptune']),
    'mass': np.array([solar_system['mercury']['mass'] / solar_system['sun']['mass'],
                      solar_system['venus']['mass'] / solar_system['sun']['mass'],
                      solar_system['earth']['mass'] / solar_system['sun']['mass'],
                      solar_system['mars']['mass'] / solar_system['sun']['mass'],
                      solar_system['jupiter']['mass'] / solar_system['sun']['mass'],
                      solar_system['saturn']['mass'] / solar_system['sun']['mass'],
                      solar_system['uranus']['mass'] / solar_system['sun']['mass'],
                      solar_system['neptune']['mass'] / solar_system['sun']['mass']]),
    'period': np.array([solar_system['mercury']['radius'] / solar_system['earth']['radius'],
                        2.4218233, 
                        4.04961,
                        6.099615,   
                        9.20669,  
                        12.35294,
                        0,
                        0,
                        0]),   # Value at index 9 is supposed to be the velocity of the celestial object
    'radius': np.array([solar_system['mercury']['radius'] / solar_system['earth']['radius'],
                        solar_system['venus']['radius'] / solar_system['earth']['radius'],
                        solar_system['earth']['radius'] / solar_system['earth']['radius'],
                        solar_system['mars']['radius'] / solar_system['earth']['radius'],
                        solar_system['jupiter']['radius'] / solar_system['earth']['radius'],
                        solar_system['saturn']['radius'] / solar_system['earth']['radius'],
                        solar_system['uranus']['radius'] / solar_system['earth']['radius'],
                        solar_system['neptune']['radius'] / solar_system['earth']['radius']]),
    'semi-major': np.array([ 11.11,  15.21,  21.44,
                             28.17,  37.1 ,  45.1 ]),
    'eccentricity': np.array([solar_system['mercury']['eccentricity'],
                              solar_system['venus']['eccentricity'],
                              solar_system['earth']['eccentricity'],
                              solar_system['mars']['eccentricity'],
                              solar_system['jupiter']['eccentricity'],
                              solar_system['saturn']['eccentricity'],
                              solar_system['uranus']['eccentricity'],
                              solar_system['neptune']['eccentricity'],])
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
