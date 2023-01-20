import numpy as np

G_meters = 6.67408e-11
G = G_meters*10**-9
J2_e = 1.082635854e-3

atm = np.array([[63.096,2.059e-4], [251.189,5.909e-11], [1000,3.561e-11]]) # atmosphereic density [altitude(km), rho]

Earth = {
    'name':'Earth',
    'mass':5.972e24,
    'mu':5.972e24*G,
    'radius':6378.0,
    'J2': 1.082635854e-3,
    'J3': -2.33936e-3*J2_e,
    'J4': -1.49601e-3*J2_e,
    'J5': -.20995e-3*J2_e,
    'J6': 0.49941e-3*J2_e,
    'zs': atm[:,0],
    'rhos': atm[:,1]*10**8, # kg/km^3
    'atm_rot_vector': np.array([0,0,72.9211e-6]), # radians per second
    'spice_file': '/home/mcmackinp/Documents/orbitalMechanics/SpiceData/de432s.bsp'
}

Sun = {
    'name':'Sun',
    'mass': 1.989e30,
    'mu':51.32712e11,
    'radius':695700.0,
    'spice_file': '/home/mcmackinp/Documents/orbitalMechanics/SpiceData/de432s.bsp'
    # 'J2': -999999999
}

Moon ={
    'name':'Moon',
    'mass':7.34767309e22,
    'mu':7.34767309e22*G, # Km^3/s^2
    'radius':1738.1,
    'J2': 202.7e-6,
    'dist2earth': 384400,
    'spice_file': '/home/mcmackinp/Documents/orbitalMechanics/SpiceData/de432s.bsp'
}