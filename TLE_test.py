import numpy as np
import matplotlib.pyplot as plt 

from OrbitPropagator import OrbitPropagaotr as OP
from OrbitPropagator import null_perts
import planetData as pData
import orbitTools as t

cb = pData.Earth

tspan = 3600*24
dt = 10

perts = null_perts()
perts['J2'] = True
perts['J3'] = True
perts['J4'] = True
perts['J5'] = True
perts['J6'] = True
perts['aero'] = True

op0 = OP(t.tle2coes('ISS_tle.txt'), tspan, dt, coes=True, deg=False, perts = perts, Cd = 2.2, A = (1e-3)**2/4)
op1 = OP(t.tle2coes('mmsats_tle.txt'), tspan, dt, coes=True, deg=False, perts = perts)
op2 = OP(t.tle2coes('progress_tle.txt'), tspan, dt, coes=True, deg=False, perts = perts)

t.plot_NOrbit(cb, [op0.rs, op1.rs, op2.rs], titles=['ISS', 'MMSATS', 'Progress'] )
