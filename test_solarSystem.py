import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
import planetData as pData
import orbitTools as tools
import spice_tools as st
from OrbitPropagator import OrbitPropagaotr as OP
from OrbitPropagator import null_perts

cb = pData.Earth

STEPS = 100000 # total steps of ephemeris data
FRAME = 'ECLIPJ2000'
OBSERVER='SUN'

# Load metadata for solar system ephemeris
spice.furnsh('SpiceData/solar_system_kernel.mk')

# get data from ephemeris file
ids,names,tcs_sec,tcs_cal=st.get_objects('SpiceData/de125.bsp', display=False)

# get only barycenters
names = [f for f in names if 'BARYCENTER' in f]

# time array for ephemeris data
times = st.tc2array(tcs_sec[0], STEPS)

# ephemeris data
rs=[]

for name in names:
    # add ephemeris data for each body to list
    rs.append(st.get_ephemeris_data(name,times,FRAME,OBSERVER))

tools.plot_NOrbit(bodies = rs, titles = names, cb = pData.Sun)
