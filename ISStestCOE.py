import numpy as np
import matplotlib.pyplot as plt 

from OrbitPropagator import OrbitPropagaotr as OP
import planetData as pData
import orbitTools as t

cb = pData.Earth

tspan = 3600*24*1
dt = 1

# ISS
c0 = [cb['radius']+414, 0.0006189, 51.6393, 0.0, 234.1955, 105.6372]

# GEO
c1 = [cb['radius']+35800, 0, 0, 0, 0, 0]

# random
c2 = [cb['radius']+3000, 0.3 ,20, 0, 15, 40]

op0 = OP(c0, tspan, dt, coes = True)
op1 = OP(c1, tspan, dt, coes = True)
op2 = OP(c2, tspan, dt, coes = True)

# op0.propagate_orbit()
# op1.propagate_orbit()
# op2.propagate_orbit()

t.plot_NOrbit(cb, [op0.rs, op1.rs,op2.rs], titles=['ISS', 'GEO', 'random'])