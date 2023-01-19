import numpy as np
import matplotlib.pyplot as plt 

from OrbitPropagator import OrbitPropagaotr as OP
import planetData as pData
import orbitTools as tool
from OrbitPropagator import null_perts

cb = pData.Earth

tspan = 3600*24
dt = 100

perts = null_perts()
perts['thrust'] = 1 # newtons
perts['thrust_direction'] = 1 # 1 in v diection, -1 opposite direction of thrust
perts['isp'] = 4300 #(s)

# initial mass 
mass0 = 50.0 #kg

# state
ra = 300+cb['radius']
rp = 215+cb['radius']

raan = 340
i = 65.2
aop = 58
ta = 332

a = (rp+ra)/2
e = (ra-rp)/(ra+rp)

state0 = [a,e,i,ta,aop,raan]

op = OP(state0, tspan, dt, coes=True, mass0=mass0, cb=cb,perts=perts)
op.plotOrbit()