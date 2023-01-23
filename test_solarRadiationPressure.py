import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
import planetData as pData
import orbitTools as tools
import spice_tools as st
from OrbitPropagator import OrbitPropagaotr as OP
from OrbitPropagator import null_perts

tspan = 3600*24*50
dt = 50

cb = pData.Earth

date0 = '2022-01-20'

# Initial conditions
iss_coes = tools.tle2coes('ISS_tle.txt')

state0 = [42095, 0.81818, 28.5, 180, 298.225, 357.857, iss_coes[-1]] # geostationary sat

perts = null_perts()

# add solar pressure perturbation
perts['srp'] = True
perts['A_srp'] = 30e-3*35e-3 #km^2
perts['CR'] = 1
mass0 = 70000

opGEO = OP(state0, tspan, dt, coes=True, perts = perts, date0 = date0, mass0 = mass0)
opISS = OP(iss_coes, tspan, dt, coes=True, perts = perts, date0 = date0, mass0 = mass0)


tools.plot_NOrbit(pData.Earth, [opISS.rs, opGEO.rs], titles = ['ISS', 'GEO'], figTitle = 'Solar radiation pressure test' )

