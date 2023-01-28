import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
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
        i*=np.pi/180
        ta*=np.pi/180
        aop*=np.pi/180
        raan*=np.pi/180
        
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
# ========= Time of Flight ==========
G_meters = 6.67408e-11
G = G_meters*10**-9
mu = 5.972e24*G

r0 = 6378+320
rm = 405081.4697

v0 = np.sqrt(2*mu*(1/r0 - 2/rm))
print(f"V_0 = {v0:.2f} km/hr")
E0 = v0**2/2 - mu/r0 # energy not eccentric anomaly
e = (rm-r0)/(rm+r0)
h = r0*v0
# tof = (2*np.pi)/mu**2 * (h/np.sqrt(1-e**2))**3/2
# a = -mu/(2*E0)
a = r0/(1-e)
tof = np.sqrt(a**3/mu) *np.pi

# Full time of flight equation
k = 0 # number of times around periapsis
E = 180 * np.pi/180 # final anomoly
Eo = 0 # starting anomoly

tof = np.sqrt(a**3/mu) * (2*np.pi*k +(E-e*np.sin(E)) - (Eo-e*np.sin(Eo)) )

print(f"Time of flight = {tof:.3}s ({(tof/3600):.3f} hours)")


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
print(getE(M[500],e))

beta = e/(1+np.sqrt(1-e**2))
nu = E + 2*np.arctan2(beta*np.sin(E),1-beta*np.cos(E))

# nu = np.arccos( ( np.cos(E) - e )/( 1-e*np.cos(E) ) ) # True anomaly

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
# ax.plot(t, E, label='E')
# ax.plot(t, nu, label='nu')
# ax.legend()

r_norm = a*(1-e**2)/(1+e*np.cos(nu))
Xpos = r_norm*np.cos(nu)
Ypos = r_norm*np.sin(nu)

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

ax.plot(r[0,:], r[1,:], r[2,:], label='S/C')
ax.plot(moonStates[:,0], moonStates[:,1], moonStates[:,2], label='Moon')
plotMoon(ax, moonStates[-1,:])
plotEarth(ax)
ax.set_aspect('equal')
ax.legend()

plt.show()