import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
import planetData as pData
import orbitTools as tools
import spice_tools as st
from OrbitPropagator import OrbitPropagaotr as OP
from OrbitPropagator import null_perts

tspan = 3600*24*25
dt = 50

cb = pData.Earth

date0 = '2022-01-20'

# Initial conditions
iss_coes = tools.tle2coes('ISS_tle.txt')

state0 = [42164, 0.001, 0, 0, 0, 0, iss_coes[-1]] # geostationary sat
# state0 = [.4055e6, 0.0549, 28.58, 0, 0, 0, iss_coes[-1]] # sat # a,e,i,ta,aop,raan

perts = null_perts()

# add n-bodies
perts['nBody'] = [pData.Moon]

opGEO = OP(state0, tspan, dt, coes=True, perts = perts, date0 = date0)
opISS = OP(iss_coes, tspan/10, dt, coes=True, perts = perts, date0 = date0)

# tools.plot_NOrbit(pData.Earth, [opISS.rs, opGEO.rs], titles = ['ISS', 'GEO'], figTitle = 'N-body w/ moon' )

tools.plot_NOrbit(pData.Earth, [opISS.rs, opGEO.rs, opGEO.perts['nBody'][0]['states'][:,:3]], titles = ['ISS', 'GEO', 'Moon'], figTitle = 'N-body w/ moon' )

