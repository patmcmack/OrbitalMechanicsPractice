import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
import planetData as pData
import orbitTools as tools
import spice_tools as st
from OrbitPropagator import OrbitPropagaotr as OP
from OrbitPropagator import null_perts

cb = pData.Earth

STEPS = 10000 # total steps of ephemeris data
FRAME = 'ECLIPJ2000'
OBSERVER='Earth'
DAYS = 22

# Load metadata for solar system ephemeris
spice.furnsh('SpiceData/solar_system_kernel.mk')

# get data from ephemeris file
ids,names,tcs_sec,tcs_cal=st.get_objects('SpiceData/de432s.bsp', display=False)

# get only barycenters
names = ['MOON']

# time array for ephemeris data
# times = st.tc2array(tcs_sec[0], STEPS)
times =np.zeros((STEPS,1))
t0 = tcs_sec[0][0] + (60*60*24*365*73.12)
t1 = tcs_sec[0][0]+(DAYS*24*60*60) + (60*60*24*365*73.12)
times[:,0] = np.linspace(t0, t1, STEPS)

# Human readable time format
tc_cal = [spice.timout(f, "YYYY MON DD HR:MN:SC") for f in [t0, t1]]

figTitle = f"From {tc_cal[0]} to {tc_cal[1]}"

# ephemeris data
rs=[]

for name in names:
    # add ephemeris data for each body to list
    rs.append(st.get_ephemeris_data(name,times,FRAME,OBSERVER))

tools.plot_NOrbit(bodies = rs, titles = names, cb = pData.Earth, figTitle = figTitle)
