import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
import math as m
from mpl_toolkits.mplot3d import Axes3D
def getE(M,e):
    # Use bisection method for stability when e = 1
    tol = 1e-8
    error = 99
    i = 1
    Ea = 0-0.1
    Eb = np.pi*2+0.1

    while i<200:
        Ec = (Ea+Eb)/2
        f_c = (Ec - e*np.sin(Ec)) - M
        f_a = (Ea - e*np.sin(Ea)) - M
        f_b = (Eb - e*np.sin(Eb)) - M
        if (abs(f_c) < tol):
            # print(f"Total iterations: {i}")
            return Ec
        elif np.sign(f_a) == np.sign(f_c):
            Ea = Ec 
        elif np.sign(f_b) == np.sign(f_c):
            Eb = Ec
        i += 1

    print("Max iterations reached")
    return -9999
        
def eci2perif(raan, aop,i):
    # Inertial to perifocal rotation matrix
    row0 = [-np.sin(raan)*np.cos(i)*np.sin(aop)+np.cos(raan)*np.cos(aop), np.cos(raan)*np.cos(i)*np.sin(aop)+np.sin(raan)*np.cos(aop), np.sin(i)*np.sin(aop)]
    row1 = [-np.sin(raan)*np.cos(i)*np.cos(aop)-np.cos(raan)*np.sin(aop), np.cos(raan)*np.cos(i)*np.cos(aop)-np.sin(raan)*np.sin(aop), np.sin(i)*np.cos(aop)]
    row2 = [np.sin(raan)*np.sin(i), -np.cos(raan)*np.sin(i), np.cos(i)]
    return np.array([row0, row1, row2])

def getMoondata(tspan):
    target = 'Moon'
    frame = 'J2000'
    observer = 'Earth'
    return np.array(spice.spkezr(target,tspan,frame,'NONE',observer)[0])

def rv2coes(state, et=0, mu=398576.0576, deg = False):
    '''
    # Get classical orbital elements from state
        Output:
    RP      Perifocal distance.
    ECC     Eccentricity.
    INC     Inclination.
    LNODE   Longitude of the ascending node.
    ARGP    Argument of periapsis.
    M0      Mean anomaly at epoch.
    T0      Epoch.
    MU      Gravitational parameter.
    NU      True anomaly at epoch.
    A       Semi-major axis. A is set to zero if
            it is not computable.
    TAU     Orbital period. Applicable only for
            elliptical orbits. Set to zero otherwise.
    '''
    
    rp,e,i,raan,aop,ma,t0,mu,ta,a,T = spice.oscltx(state,et,mu)

    if deg:
        i/=np.pi/180
        ta/=np.pi/180
        aop/=np.pi/180
        raan/=np.pi/180
        
    return[a,e,i,ta,aop,raan]

def plotEarth(ax):
    radius = 6378.0

    u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:30j]
    x = radius*np.cos(u)*np.sin(v)
    y = radius*np.sin(u)*np.sin(v)
    z = radius*np.cos(v)
    ax.plot_surface(x,y,z,cmap = 'Blues', zorder = 0.3, alpha = 0.1, edgecolors=[.5,.5,.5], linewidth=0.1)

def plotMoon(ax, moonstate):
    radius = 1738.1

    u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:30j]
    x = radius*np.cos(u)*np.sin(v) + moonstate[0]
    y = radius*np.sin(u)*np.sin(v) + moonstate[1]
    z = radius*np.cos(v) + moonstate[2]
    ax.plot_surface(x,y,z,cmap = 'Greys', zorder = 0.3, alpha = 0.1, edgecolors=[.5,.5,.5], linewidth=0.1)

def getPosTLI(t, e, a, Ef, raan, aop, i, mu):
    '''
    Get S/C position after TLI
    Assumed to be Hohmann transfer
    '''
    # E = Ef
    # beta = e/(1+np.sqrt(1-e**2))
    # nu =  E + 2*np.arctan2(beta*np.sin(E),1-beta*np.cos(E)) 
    # r_norm = a*(1-e**2)/(1+e*np.cos(nu))
    # # get r and v vector in perifocal frame
    # r_perif = r_norm*np.array([np.cos(nu), np.sin(nu), 0])
    # v_perif = np.sqrt(mu*a)/r_norm*np.array([-np.sin(E), np.cos(E)*np.sqrt(1-e**2),0])

    # # Rotation matrix from perifocal to ECI
    # perif2eci = np.transpose(eci2perif(raan, aop, i))

    # # calc r and v in inertial frame
    # r_sat = np.dot(perif2eci, r_perif)
    # v_sat = np.dot(perif2eci, v_perif)

    # ---------------------
    n = np.sqrt(mu/a**3) # mean motion (\dot{M})
    M = n * t
    E = np.zeros(M.shape)
    for i, m in enumerate(M):
        E[i] = getE(m,e) # Kepler's eqn

    beta = e/(1+np.sqrt(1-e**2))
    nu = E + 2*np.arctan2(beta*np.sin(E),1-beta*np.cos(E))
    r_norm = a*(1-e**2)/(1+e*np.cos(nu))

    # get r and v vector in perifocal frame
    r_perif = r_norm*np.array([np.cos(nu), np.sin(nu), np.zeros(len(nu))])
    v_perif = np.sqrt(mu*a)/r_norm*np.array([-np.sin(E), np.cos(E)*np.sqrt(1-e**2),np.zeros(len(nu))])

    # Rotation matrix from perifocal to ECI
    perif2eci = np.transpose(eci2perif(raan, aop, i))

    # calc r and v in inertial frame
    r_sat = np.dot(perif2eci, r_perif)
    v_sat = np.dot(perif2eci, v_perif)

    # return r_perif
    return r_sat

# ========= Initial Conditions =======

# Central body = Earth
G_meters = 6.67408e-11
G = G_meters*10**-9
mu = 5.972e24*G

# Time states
date0 = '2022-01-20'
dt = 100 # second
tspan = 60*60*24*100 # seconds

# Get moon data from JPL SPICE data 
print("Fetching Moon SPICE data")
spice.furnsh('SpiceData/solar_system_kernel.mk')
start_time = spice.utc2et(date0) # convert start date to seconds after J2000
t = np.linspace(start_time, start_time+tspan, m.floor(tspan/dt))
moonStates = getMoondata(t)
avgMoonRadius = np.mean(np.linalg.norm(moonStates[:,:3], axis=1))
coes_start = rv2coes(moonStates[0])
coes_end = rv2coes(moonStates[-1])
avgMoonInc = (coes_end[2]+coes_start[2])/2.0

iterations = len(t)
# Parking orbit of spacecraft (assume near lunar plane)
ra = 6378+500
rp = 6378+500
aop = 96*np.pi/180#(coes_end[4]+coes_start[4])/2.0 # argument of periapsis
i = avgMoonInc # Inclination
raan = (coes_end[5]+coes_start[5])/2.0 # right of assending node
e = (ra-rp)/(ra+rp)
a = (ra+rp)/2
n = np.sqrt(mu/a**3) # mean motion (\dot{M})
T = 2*np.pi*np.sqrt(a**3/mu) # Orbital period 

# correct for multiple orbits 
M = n * np.linspace(0, T, m.floor(T/dt)) # Mean anomaly
# M = n * ((t-start_time) - T*np.floor((t-start_time)/T))# Mean anomaly
E = np.zeros(M.shape) # Eccentric anomaly
print("Calculating S/C parking orbit in time")
for j, m in enumerate(M):
    print('\r', end='')
    print(f"\t{(j/len(M)*100):.1f} %", end='')
    E[j] = getE(m,e) # solve Kepler's Eqn
print('')

# Repeat E to get full t span
E = np.tile(E, int(np.floor(len(t)/len(E))))
E = np.concatenate((E, E[:int(len(t)-len(E))]))

# Get Spacecraft position in time in parking orbit 
beta = e/(1+np.sqrt(1-e**2))
nu = E + 2*np.arctan2(beta*np.sin(E),1-beta*np.cos(E)) # True anomoly
r_norm = a*(1-e**2)/(1+e*np.cos(nu))
# get r and v vector in perifocal frame
r_perif = r_norm*np.array([np.cos(nu), np.sin(nu), np.zeros(len(nu))])
v_perif = np.sqrt(mu*a)/r_norm*np.array([-np.sin(E), np.cos(E)*np.sqrt(1-e**2),np.zeros(len(nu))])

# Rotation matrix from perifocal to ECI
perif2eci = np.transpose(eci2perif(raan, aop, i))

# calc r and v in inertial frame
r_sat = np.dot(perif2eci, r_perif)
v_sat = np.dot(perif2eci, v_perif)

# ====== Iterate in time ======
fig = plt.figure()
ax = fig.add_subplot(111, projection ='3d')
ax.plot(r_sat[0,:], r_sat[1,:], r_sat[2,:], label='S/C parking')
# ax = Axes3D(fig)
r_min = np.empty(t.shape)
r_min[:] = np.nan
print("Calculating TLI injections")
for j, time in enumerate(t):
    # j=650
    # Iterate through each point in time 

    # Estimate time of flight (Assume Hohmann transfer)
    rp = r_norm[j] # initial altitude (rp)
    ra = np.linalg.norm(moonStates[j, :3])
    e = (ra-rp)/(ra+rp)
    a = (ra+rp)/2
    # v0 = np.sqrt(2*mu*(1/r0 - 2/avgMoonRadius)) # not used

    k = 0 # number of times around periapsis
    Ef = 180 * np.pi/180 # final anomoly
    Eo = 0 # starting anomoly

    tof = np.sqrt(a**3/mu) * (2*np.pi*k +(Ef-e*np.sin(Ef)) - (Eo-e*np.sin(Eo)) )

    # Loop end condition
    if time+tof >= t[-1]:
        break

    # ========= Moon postion after TLI =========
    idx = np.searchsorted(t, time+tof, side ='left')
    r_moon = moonStates[idx, :3]

    coes_moon = rv2coes(moonStates[idx])
    i_moon = coes_moon[2]
    raan_moon = coes_moon[5]

    # ========= Spacecraft position after TLI =========
    sat_state = np.concatenate((r_sat[:, j], v_sat[:, j])) 
    coes_sat = rv2coes(sat_state)
    nu_sat = coes_sat[3]
    aop_sat = nu_sat#coes_sat[4]
    raan_sat = coes_sat[5]
    i_sat = coes_sat[2]
    Ef = 180 *np.pi/180
    aop_j = aop_sat#96 * np.pi/180#nu[j] # is this correct???
    raan_j = raan_sat#90*np.pi/180
    t_j = t[j:idx] - t[j]
    # r_sat_TLI = getPosTLI(t_j, e, a, Ef, raan_moon, aop_j, i_moon, mu)
    # v_perif_tli = np.sqrt(mu*a)/r_norm*np.array([-np.sin(E), np.cos(E)*np.sqrt(1-e**2),np.zeros(len(nu))])
    # ------------------------------------
    n = np.sqrt(mu/a**3) # mean motion (\dot{M})
    M = n * t_j
    E = np.zeros(M.shape)
    for jj, m in enumerate(M):
        E[jj] = getE(m,e)

    beta = e/(1+np.sqrt(1-e**2))
    nu = E + 2*np.arctan2(beta*np.sin(E),1-beta*np.cos(E))
    r_norm = a*(1-e**2)/(1+e*np.cos(nu))

    r_perif = r_norm*np.array([np.cos(nu), np.sin(nu), np.zeros(len(nu))])
    v_perif = np.sqrt(mu*a)/r_norm*np.array([-np.sin(E), np.cos(E)*np.sqrt(1-e**2),np.zeros(len(nu))])

    aop = nu_sat+aop_sat#96 * np.pi/180 # argument of periapsis
    # i = 28.5 * np.pi/180

    # Rotation matrix from perifocal to ECI
    perif2eci = np.transpose(eci2perif(raan, aop, i))

    # calc r and v in inertial frame
    r_sat_TLI = np.dot(perif2eci, r_perif)
    # ------------------------------------

    

    # ========= Calc distance from moon post TLI =========
    r_moon2sc = np.linalg.norm(r_sat_TLI[:,-1] - r_moon)
    r_min[j] = r_moon2sc

    # ========= Output progress =========
    if j%(len(t)/100)<=1:
        print('\r', end='')
        print(f"\t{(j/len(t) * 100):.0f} %", end='')

    # # ax.clear()
    # ax.plot(r_sat[0,j], r_sat[1,j], r_sat[2,j], marker='o')
    # ax.plot(r_sat_TLI[0,:], r_sat_TLI[1,:], r_sat_TLI[2,:], label='S/C')
    # ax.plot(moonStates[j:idx,0], moonStates[j:idx,1], moonStates[j:idx,2], label='Moon')
    # plotMoon(ax, moonStates[j,:])
    # plotMoon(ax, moonStates[idx,:])
    # plotEarth(ax)
    # ax.set_aspect('equal')
    # ax.legend()
    # ax.set_title(j)

    # plt.show()
    print(j)
print('')

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(t, r_min)
plt.show()


# ========= Position vs time (Kepler's Problem) ==========
ra = rm
rp = r0 
t = np.linspace(0,tof,5000) # time vector
e = (ra-rp)/(ra+rp)

E = 0 # Eccentric anomaly
nu = 0 # true anomaly
n = np.sqrt(mu/a**3) # mean motion (\dot{M})
M = n * t
E = np.zeros(M.shape)
for i, m in enumerate(M):
    E[i] = getE(m,e)
# M = E - e*np.sin(E) # Mean anomaly

beta = e/(1+np.sqrt(1-e**2))
nu = E + 2*np.arctan2(beta*np.sin(E),1-beta*np.cos(E))

# nu = np.arccos( ( np.cos(E) - e )/( 1-e*np.cos(E) ) ) # True anomaly

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
# ax.plot(t, E, label='E')
# ax.plot(t, nu, label='nu')
# ax.legend()

r_norm = a*(1-e**2)/(1+e*np.cos(nu))
# Xpos = r_norm*np.cos(nu)
# Ypos = r_norm*np.sin(nu)

# ax.plot(Xpos, Ypos)


# plt.show()
# ======= 3D form =======
spice.furnsh('SpiceData/solar_system_kernel.mk')
moonStates = getMoondata(t)

# get classical orbital elements at end of flight
coes_end = rv2coes(moonStates[-1])
i = coes_end[2]
raan = coes_end[5]
# get r and v vector in perifocal frame
r_perif = r_norm*np.array([np.cos(nu), np.sin(nu), np.zeros(len(nu))])
v_perif = np.sqrt(mu*a)/r_norm*np.array([-np.sin(E), np.cos(E)*np.sqrt(1-e**2),np.zeros(len(nu))])

aop = 96 * np.pi/180 # argument of periapsis
# i = 28.5 * np.pi/180

# Rotation matrix from perifocal to ECI
perif2eci = np.transpose(eci2perif(raan, aop, i))

# calc r and v in inertial frame
r = np.dot(perif2eci, r_perif)
v = np.dot(perif2eci, v_perif)

# Calc min distance from moon
r_moon2sc = np.linalg.norm(np.transpose(r)-moonStates[:,:3], axis=1)

print(f"Min distance from moon after TLI: {np.min(r_moon2sc):.2f} km")

ax.plot(r[0,:], r[1,:], r[2,:], label='S/C')
ax.plot(moonStates[:,0], moonStates[:,1], moonStates[:,2], label='Moon')
plotMoon(ax, moonStates[-1,:])
plotEarth(ax)
ax.set_aspect('equal')
ax.legend()

plt.show()