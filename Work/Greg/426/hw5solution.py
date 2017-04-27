# PHY 494 Homework 06 --- Orbit of Uranus: Solution
# Copyright (c) 2016 Oliver Beckstein
# All Rights Reserved.
#
# This solution may not be distributed.
#
# Students of the PHY494 class may use code from this solution
# for other projects in this class.

import numpy as np


# gravitational constant in astronomical units
G_gravity = 4*np.pi**2

#============================================================
# Parameters of the problem
#------------------------------------------------------------
#
# store parameters in dict; use as e.g.
#  theta0['Uranus']
# to access initial angle of Uranus.

# angular position in 1690 in degrees
theta0 = {'Uranus': 205.640,
         'Neptune': 288.380,
}
# distance from the sun in AU
distance = {'Uranus': 19.1914,
            'Neptune': 30.0611,
}
# mass in AU
mass = {'Sun': 1.,
        'Uranus': 4.366244e-5,
        'Neptune': 5.151389e-5,
}
# orbital period in Earth years
period = {'Uranus': 84.0110,
          'Neptune': 164.7901,
}

#------------------------------------------------------------


def initial_position(angle, distance):
    """Calculate initial planet position.

    Parameters
    ----------
    angle : float
       initial angle relative to x axis (in degrees)
    distance : float
       initial distane from sun (in AU)

    Returns
    -------
    array
       position (x, y)
    """
    x = np.deg2rad(angle)
    return distance * np.array([np.cos(x), np.sin(x)])

def initial_velocity(angle, distance, period):
    x = np.deg2rad(angle)
    v = 2*np.pi*distance / period
    return v * np.array([np.sin(x), -np.cos(x)])

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
    return -G_gravity*m*M/rr * rhat

def F_planets(positions, coupled=True):
    rU = positions[0]
    rN = positions[1]
    if coupled:
        F_UN = F_gravity(rU - rN, mass['Uranus'], mass['Neptune'])
    else:
        F_UN = 0
    F_NU = -F_UN
    F_U = F_gravity(rU, mass['Uranus'], mass['Sun']) + F_UN
    F_N = F_gravity(rN, mass['Neptune'], mass['Sun']) + F_NU
    return np.array([F_U, F_N])

def omega(v, r):
    """Calculate angular velocity.

    The angular velocity is calculated as
    .. math::

          \omega = \frac{|\vec{v}|}{|\vec{r}|}

    Parameters
    ----------
    v : array
       velocity vectors for all N time steps; this
       should be a (N, dim) array
    r : array
       position vectors (N, dim) array

    Returns
    -------
    array
       angular velocity for each time step as 1D array of
       length N
    """
    speed = np.linalg.norm(v, axis=1)
    distance = np.linalg.norm(r, axis=1)
    return speed/distance

def integrate_orbits(dt=0.1, t_max=320, coupled=True):
    """Integrate equations of motions of Uranus and Neptune.

    Parameters
    ----------
    dt : float
       integrator timestep
    t_max : float
       integrate to t_max years
    coupled : bool
       * `True`: include the interaction between Neptune and Uranus
       * `False`: no interaction (Uranus and Neptune move independently)

    Returns
    -------
    time : array
       array with the times for each step (in years)
    r : array
       positions of the planets, shape (N, 2, 2), where::
          r[:, 0] = rU : Uranus x, y
          r[:, 1] = rN : Neptune x, y
    v : array
       velocities of the planets, shape (N, 2, 2), where::
          v[:, 0] = vU : Uranus vx, vy
          v[:, 1] - vN : Neptune vx, vy
    """
    nsteps = int(t_max/dt)
    time = dt * np.arange(nsteps)

    # shape = (step, planet, x/y)
    r = np.zeros((nsteps, 2, 2))
    v = np.zeros_like(r)

    r[0, 0, :] = initial_position(theta0['Uranus'], distance['Uranus'])
    r[0, 1, :] = initial_position(theta0['Neptune'], distance['Neptune'])
    v[0, 0, :] = initial_velocity(theta0['Uranus'],
                                  distance['Uranus'], period['Uranus'])
    v[0, 1, :] = initial_velocity(theta0['Neptune'],
                                  distance['Neptune'], period['Neptune'])
    # note: masses has to be shape (2,1) for broadcasting
    #       for F/masses
    masses = np.array([[mass['Uranus']], [mass['Neptune']]])

    # start force evaluation for first step
    Ft = F_planets(r[0], coupled=coupled)
    for i in range(nsteps-1):
        vhalf = v[i] + 0.5*dt * Ft/masses
        r[i+1, :] = r[i] + dt * vhalf
        Ftdt = F_planets(r[i+1], coupled=coupled)
        v[i+1, :] = vhalf + 0.5*dt * Ftdt/masses
        # new force becomes old force
        Ft = Ftdt

    return time, r, v


if __name__ == "__main__":
    # The following code only runs when the script is executed from the shell with
    #
    #    python ./outerplanets.py
    #
    # or inside a notebook with
    #
    #    %run ./outerplanets.py
    #
    # If you 'import outerplanets' then only the function definitions will be
    # executed, and *not* the following code.
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.style.use("ggplot")


    print("Simulating Uranus and Neptune orbits WITHOUT interactions")
    time, r0, v0 = integrate_orbits(t_max=30, coupled=False)
    rU0 = r0[:, 0]
    vU0 = v0[:, 0]
    omegaU0 = omega(vU0, rU0)

    print("Simulating Uranus and Neptune orbits WITH interactions")
    time, r, v = integrate_orbits(t_max=30, coupled=True)
    rU = r[:, 0]
    rN = r[:, 1]
    vU = v[:, 0]
    omegaU = omega(vU, rU)

    #----------
    # IMPLEMENT
    #
    DeltaOmegaU = omegaU0 - omegaU
    #
    #---------

    # plot orbits
    #fig_orbits = "uranus_neptune_orbits.pdf"
    #fig_anomaly = "uranus_anomaly.pdf"

    ax = plt.subplot()
    ax.plot(rU[:,0], rU[:, 1], label="Uranus")
    ax.plot(rN[:,0], rN[:, 1], label="Neptune")
    ax.plot(rU0[:,0], rU0[:,1], linestyle="--", label="Uranus (no Neptune)")
    ax.set_aspect(1)
    ax.set_xlabel(r"$x$ (AU)")
    ax.set_ylabel(r"$y$ (AU)")
    ax.legend(loc="upper right")
    ax.set_title("Uranus and Neptune orbits")
    plt.show()
""" ax = plt.subplot(1,2,2)
    ax.plot(time, DeltaOmegaU)
    ax.set_xlabel("years")
    ax.set_ylabel(r"Uranus anomaly $\Delta\omega_U$")
    ax.set_title("Uranus anomaly")
    ax.figure.tight_layout()
    ax.figure.savefig(fig_anomaly)
    print("Uranus anomaly plotted in {0}".format(fig_anomaly))

"""
print(r)
print(r)[0]
print(r)[0,0]