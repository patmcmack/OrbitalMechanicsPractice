import numpy as np
import matplotlib.pyplot as plt 

from OrbitPropagator import OrbitPropagaotr as OP
import planetData as pData
import orbitTools as t

cb = pData.Earth

# Initial conditions 
r_mag = cb['radius'] + 2500
v_mag = np.sqrt(cb['mu']/r_mag)

r0 = [r_mag, 0, 0]
v0 = [v_mag*.01, v_mag*1.25 , v_mag*.15]

r02 = [r_mag, 0, 0]
v02 = [v_mag*.01, v_mag*1 , v_mag*.1]

tspan = 3600*24*2

dt = 10

op = OP(r0,v0,tspan,dt)
op.propagate_orbit()

# op.plotOrbit()

op2 = OP(r02,v02,tspan,dt)
op2.propagate_orbit()

# t.plot_NOrbit(cb, [op.rs], titles=['1'])
t.plot_NOrbit(cb, [op.rs,op2.rs], titles=['1','2'])
