import numpy as np

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
            print(f"Total iterations: {i}")
            return Ec
        elif np.sign(f_a) == np.sign(f_c):
            Ea = Ec 
        elif np.sign(f_b) == np.sign(f_c):
            Eb = Ec
        i += 1

    print("Max iterations reached")
    return -9999
        


# ========= Time of Flight ==========
G_meters = 6.67408e-11
G = G_meters*10**-9
mu = 5.972e24*G

r0 = 6378+320
rm = 384400

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
t = np.linspace(0,tof,1000) # time vector
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