import numpy as np
import matplotlib.pyplot as plt
import math as m
import datetime

import planetData as pData
import spiceypy as spice #https://spiceypy.readthedocs.io/en/v2.3.1/documentation.html

def plot_NOrbit(cb, bodies, titles, figTitle = ''):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')

        u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:30j]
        x = cb['radius']*np.cos(u)*np.sin(v)
        y = cb['radius']*np.sin(u)*np.sin(v)
        z = cb['radius']*np.cos(v)
        ax.plot_surface(x,y,z,cmap = 'Blues', zorder = 0.3, alpha = 0.1, edgecolors=[.5,.5,.5], linewidth=0.1)

        maxVal = 0
        for n, b in enumerate(bodies):
            ax.plot(b[:,0], b[:,1], b[:,2], linestyle='-',label=titles[n],  zorder = 0.5)
            ax.plot(b[0,0], b[0,1], b[0,2], marker='o')
            if np.max(np.abs(b)) > maxVal:
                maxVal = np.max(np.abs(b))

        

        # ax.set_xlim([-maxVal, maxVal])
        # ax.set_ylim([-maxVal, maxVal])
        # ax.set_zlim([-maxVal, maxVal])

        ax.set_xlabel(['X (km)'])
        ax.set_ylabel(['Y (km)'])
        ax.set_zlabel(['Z (km)'])
        if len(figTitle)>0:
            ax.set_title(figTitle)

        ax.set_aspect('equal')
        ax.legend()

        # ax.set_title(title)

        plt.show()

def coes2rv(coes, deg=False, mu = pData.Earth['mu']):
    # Convert Classical orbital elements to state (r & v)
    
    a,e,i,ta,aop,raan, date = coes
    # a,e,i,ta,aop,raan = coes
    if deg:
        i*=np.pi/180
        ta*=np.pi/180
        aop*=np.pi/180
        raan*=np.pi/180
    E = ecc_anomaly([ta,e], 'tae')
    r_norm=a*(1-e**2)/(1+e*(np.cos(ta)))

    # get r and v vector in perifocal frame
    r_perif = r_norm*np.array([m.cos(ta), m.sin(ta),0])
    v_perif = m.sqrt(mu*a)/r_norm*np.array([-m.sin(E), m.cos(E)*m.sqrt(1-e**2),0])

    # Rotation matrix from perifocal to ECI
    perif2eci = np.transpose(eci2perif(raan, aop, i))

    # calc r and v in inertial frame
    r = np.dot(perif2eci, r_perif)
    v = np.dot(perif2eci, v_perif)

    return r,v 

def eci2perif(raan, aop,i):
    # Inertial to perifocal rotation matrix
    row0 = [-m.sin(raan)*m.cos(i)*m.sin(aop)+m.cos(raan)*m.cos(aop), m.cos(raan)*m.cos(i)*m.sin(aop)+m.sin(raan)*m.cos(aop), m.sin(i)*m.sin(aop)]
    row1 = [-m.sin(raan)*m.cos(i)*m.cos(aop)-m.cos(raan)*m.sin(aop), m.cos(raan)*m.cos(i)*m.cos(aop)-m.sin(raan)*m.sin(aop), m.sin(i)*m.cos(aop)]
    row2 = [m.sin(raan)*m.sin(i), -m.cos(raan)*m.sin(i), m.cos(i)]
    return np.array([row0, row1, row2])

def getDensity(alt):
    # Estimate density at a specific altitude in Km
    rhos, zs = find_rho_z(alt)
    if rhos[0]==0:
        return 0
    Hi = -(zs[1]-zs[0])/m.log(rhos[1]/rhos[0])

    return rhos[0]*m.exp(-(alt-zs[0])/Hi)

def find_rho_z(z, zs = pData.Earth['zs'], rhos = pData.Earth['rhos']):
    if not 1<z<1000:
        return [[0,0],[0,0]]

    # Get two surrounding data points 
    for n in range(len(rhos)-1):
        if zs[n]<z<zs[n+1]:
            return [[rhos[n], rhos[n+1]], [zs[n], zs[n+1]]]
    
    # out of range
    return [[0,0],[0,0]]

def tle2coes(fileName, mu = pData.Earth['mu']):
    # read TLE file and return classical orbital elements

    with open(fileName,'r') as f:
        lines = f.readlines()
    
    line0 = lines[0].strip()
    line1 = lines[1].strip().split()
    line2 = lines[2].strip().split()

    epoch = line1[3]
    year, month, day, hour = calc_epoch(epoch)

    # COE
    i = float(line2[2])*np.pi/180
    raan = float(line2[3])*np.pi/180
    e = line2[4]
    e = float('0.'+e)
    aop = float(line2[5])*np.pi/180
    Me = float(line2[6])*np.pi/180 # radians
    mean_motion = float(line2[7]) #revs/day
    T = 1/mean_motion*24*3600 # seconds (Period)
    a = (T**2*mu/4/np.pi**2)**(1/3)

    E = ecc_anomaly([Me,e],'newton')
    ta = true_anomaly([E,e])

    return a,e,i,ta,aop,raan,[year,month,day,hour]


def calc_epoch(epoch):
    year = int('20'+epoch[:2])

    epoch = epoch[2:].split('.')

    # day of year
    doy = int(epoch[0])-1

    hour = float('0.'+epoch[1])*24

    date = datetime.date(year,1,1)+datetime.timedelta(doy)

    month = float(date.month)
    day = float(date.day)

    return year, month, day, hour

def ecc_anomaly(arr,method,tol=1e-8):
    # Calculate eccentric anomaly (E)
    if method =='newton':
        Me,e = arr
        if Me<np.pi/2.0: 
            E0=Me+e/2.0
        else:
            E0 = Me-e 
        for n in range(200): # Max num of iterations
            ratio = (E0-e*np.sin(E0)-Me)/(1-e*np.cos(E0))
            if abs(ratio)<tol:
                if n==0:
                    return E0
                else:
                    return E1 
            else:
                E1 = E0-ratio
                E0=E1
        return False # did not converg
    elif method=='tae':
        ta, e = arr
        E = 2*m.atan(m.sqrt((1-e)/(1+e))*m.tan(ta/2.0))
        return E
    else:
        print('Invalid method passed to ecc_anomaly()')
        return False

def coes2rvSPICE(coes, et=0, T0=0, deg=False, mu = pData.Earth['mu']):
    '''
    Spice implementation of classical orbital elements to state
        INPUT:
    RP      Perifocal distance.
    ECC     Eccentricity.
    INC     Inclination.
    LNODE   Longitude of the ascending node.
    ARGP    Argument of periapse.
    M0      Mean anomaly at epoch.
    T0      Epoch.
    MU      Gravitational parameter.
    '''
    rp,e,i,raan,aop,ma = coes

    if deg:
        i*=np.pi/180
        ta*=np.pi/180
        aop*=np.pi/180
        raan*=np.pi/180
    state = spice.conics([rp,e,i,raan,aop,ma, T0, mu], et)
    return(state)

def rv2coes(state, et=0, mu=pData.Earth['mu'], deg = True):
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

def true_anomaly(arr):
    E,e = arr
    return 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2.0))

'''
spice.utc2et(self.date0) # seconds after J2000
'''