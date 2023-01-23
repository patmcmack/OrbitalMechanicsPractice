import matplotlib.pyplot as plt
from scipy.integrate import ode
import numpy as np
import orbitTools as tools
import spiceypy as spice
import spice_tools as st

import planetData as pData

def null_perts(): # Null Perturbations
    return {
        'J2': False,
        'J3': False,
        'J4': False,
        'J5': False,
        'J6': False,
        'aero': False,
        'thrust': 0,
        'thrust_direction': 0,
        'isp': 0,
        'nBody': [],
        'srp': False,
        'CR': 0, # coefficient of reflectivity
        'A_srp': 0 # area for solar radiation pressure 
    }

class OrbitPropagaotr:
    def __init__(self, state0, tspan, dt, coes = False, mass0 = 99999, deg=True, cb = pData.Earth, perts = null_perts(), Cd = 0, A = 0, date0 = '2022-01-01', frame='J2000'):
        if coes:
            self.r0,self.v0 = tools.coes2rv(state0, deg = deg, mu = cb['mu'])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]

        
        self.tspan = tspan
        self.dt = dt 
        self.cb = cb
        self.step = 1
        self.perts = perts
        self.Cd = Cd # coef of drag
        self.A = A # cross-sectional area
        self.mass0 = mass0
        self.date0 = date0
        self.frame = frame

        self.nsteps = int(np.ceil(self.tspan/self.dt))

        # Initial state
        self.y0 = self.r0.tolist() + self.v0.tolist()+[self.mass0]
        self.ts = np.zeros((self.nsteps,1))
        self.ys = np.zeros((self.nsteps,7))
        self.ts[0] = 0
        self.ys[0,:] = self.y0 
        self.spice_files_loaded = []

        # Check if loading in spice data
        if self.perts['nBody'] or self.perts['srp']:
            # load kernel and log
            spice.furnsh('SpiceData/solar_system_kernel.mk')
            self.spice_files_loaded.append('SpiceData/solar_system_kernel.mk')

            # convert start date to seconds after J2000
            self.start_time = spice.utc2et(self.date0)

            self.spice_tspan = np.linspace(self.start_time, self.start_time+self.tspan, self.nsteps)

            # get sun state if srp
            if self.perts['srp']:

                spice.furnsh(self.cb['spice_file'])

                self.spice_files_loaded.append(self.cb['spice_file'])

                # calc central body states throught propagation with respect to sun
                self.cb['states']=st.get_ephemeris_data(self.cb['name'], self.spice_tspan, self.frame, 'SUN')

        # Load kernels for each body
        for body in self.perts['nBody']:

            if body['spice_file'] not in self.spice_files_loaded:
                spice.furnsh(body['spice_file'])
                self.spice_files_loaded.append(body['spice_file'])

            # Calculate body states with respect to central body
            body['states'] = st.get_ephemeris_data(body['name'],self.spice_tspan,self.frame, self.cb['name'])

        # solver 
        self.solver = ode(self.diffy_q)
        # self.solver.set_integrator('lsoda') # Adams Method/BDF
        self.solver.set_integrator('dopri5') # Runge-Kutta (4)5 method
        self.solver.set_initial_value(self.y0,0)

        self.propagate_orbit()

    def propagate_orbit(self):
        
        # step through orbit
        print("Propagating orbit:")
        while self.solver.successful() and self.step < self.nsteps:
            

            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step] = self.solver.t
            self.ys[self.step] = self.solver.y
            self.step += 1
            if self.step%(self.nsteps/100)<=1:
                print('\r', end='')
                print(f"\t{(self.step/self.nsteps * 100):.1f} %", end='')
        print('')

        self.rs = self.ys[:,:3]
        self.vs = self.ys[:,3:]

    def diffy_q(self, t, y):
        rx,ry,rz,vx,vy,vz,mass = y
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])
        dmdt = 0

        norm_r = np.linalg.norm(r)

        a = -r*self.cb['mu']/norm_r**3

        # J2 perturbations
        if self.perts['J2']:
            # Orbital Mechanics of Eng. Students H. Curtis 4th ed. Eqn. 10.30
            z2 = r[2]**2
            r2 = norm_r**2
            tx = r[0]/norm_r*(5*z2/r2-1)
            ty = r[1]/norm_r*(5*z2/r2-1)
            tz = r[2]/norm_r*(5*z2/r2-3)

            a_j2 = 1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])

            a+=a_j2

        # J3 perturbations
        if self.perts['J3']:
            #Schaub, H., Junkins, J.L., 2009. Analytical Mechanics of Space Systems, second ed (page 553)
            tx = 5*(7*(r[2]/norm_r)**3 - 3*(r[2]/norm_r))
            ty = 5*(7*(r[2]/norm_r)**3 - 3*(r[2]/norm_r))
            tz = 3*(1-10*(r[2]/norm_r)**2 + 35/3*(r[2]/norm_r)**4)

            a_j3 = 0.5*self.cb['J3']*(self.cb['mu']/norm_r**2)*(self.cb['radius']/norm_r)**3*np.array([tx,ty,tz])

            a+=a_j3

        # J4 perturbations
        if self.perts['J4']:
            #Schaub, H., Junkins, J.L., 2009. Analytical Mechanics of Space Systems, second ed (page 553)
            tx = (3-42*(r[2]/norm_r)**2 + 63*(r[2]/norm_r)**4)
            ty = (3-42*(r[2]/norm_r)**2 + 63*(r[2]/norm_r)**4)
            tz = (15-70*(r[2]/norm_r)**2 + 63*(r[2]/norm_r)**4)

            a_j4 = (5/8)*self.cb['J4']*(self.cb['mu']/norm_r**2)*(self.cb['radius']/norm_r)**4*np.array([tx,ty,tz])

            a+=a_j4

        # J5 perturbations
        if self.perts['J5']:
            #Schaub, H., Junkins, J.L., 2009. Analytical Mechanics of Space Systems, second ed (page 553)
            tx = 3*(35*(r[2]/norm_r) - 210*(r[2]/norm_r)**3 + 231*(r[2]/norm_r)**5)
            ty = 3*(35*(r[2]/norm_r) - 210*(r[2]/norm_r)**3 + 231*(r[2]/norm_r)**5)
            tz = (693*(r[2]/norm_r)**6 - 945*(r[2]/norm_r)**4 + 315*(r[2]/norm_r)**2 - 15)

            a_j5 = (1/8)*self.cb['J5']*(self.cb['mu']/norm_r**2)*(self.cb['radius']/norm_r)**5*np.array([tx,ty,tz])

            a+=a_j5

        # J6 perturbations
        if self.perts['J6']:
            #Schaub, H., Junkins, J.L., 2009. Analytical Mechanics of Space Systems, second ed (page 553)
            tx = (35 - 945*(r[2]/norm_r)**2 + 3465*(r[2]/norm_r)**4 - 3003*(r[2]/norm_r)**6)
            ty = (35 - 945*(r[2]/norm_r)**2 + 3465*(r[2]/norm_r)**4 - 3003*(r[2]/norm_r)**6)
            tz = (245 - 2205*(r[2]/norm_r)**2 + 4851*(r[2]/norm_r)**4 - 3003*(r[2]/norm_r)**6)

            a_j6 = (-1/16)*self.cb['J6']*(self.cb['mu']/norm_r**2)*(self.cb['radius']/norm_r)**6*np.array([tx,ty,tz])

            a+=a_j6

        # Aerodynamic drag perturbation
        if self.perts['aero']:
            alt = norm_r - self.cb['radius']
            rho = tools.getDensity(alt)

            # velocity relative to rotating atmosphere
            v_rel = v-np.cross(self.cb['atm_rot_vector'],r)

            a_drag = -0.5 * rho * v_rel * np.linalg.norm(v_rel)*self.Cd*self.A/mass

            a+=a_drag

        # Thrust perturbation
        if self.perts['thrust']:
            # thrust vector
            a+=self.perts['thrust_direction']*(np.array(v)/np.linalg.norm(v))*self.perts['thrust']/mass/1000 # km/s^2

            # Derivative of total mass
            dmdt =- self.perts['thrust']/self.perts['isp']/9.81

        # N-body perturbation
        for body in self.perts['nBody']:
            # Third-Body Perturbation Effects on Satellite Formations, C. Roscoe, 2015

            # central body to n-body
            r_cb2nb = body['states'][self.step,:3]

            # satellite to body
            r_sat2body = r_cb2nb-r

            a +=  body['mu'] * (r_sat2body/np.linalg.norm(r_sat2body)**3 - r_cb2nb/(np.linalg.norm(r_cb2nb))**3)

        # Solar radiation pressure
        if self.perts['srp']:
            # ORBIT MECHANICS ABOUT SMALL ASTEROIDS, D.J. Scheeres, 2007

            # vector form sun to spacecraft 
            r_sun2sc = self.cb['states'][self.step, :3]+r 

            a += (1+self.perts['CR'])*pData.Sun['G1']*self.perts['A_srp']/mass/np.linalg.norm(r_sun2sc)**3*r_sun2sc

        return[vx,vy,vz,a[0],a[1],a[2], dmdt]

    def plotOrbit(self, title='Title'):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')

        u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:30j]
        x = self.cb['radius']*np.cos(u)*np.sin(v)
        y = self.cb['radius']*np.sin(u)*np.sin(v)
        z = self.cb['radius']*np.cos(v)
        ax.plot_surface(x,y,z,cmap = 'Blues', zorder = 0.3, alpha = 0.1, edgecolors=[.5,.5,.5], linewidth=0.1)

        ax.plot(self.rs[:,0], self.rs[:,1], self.rs[:,2], 'k', label='Tragectory', zorder = 0.5)
        ax.plot(self.rs[0,0], self.rs[0,1], self.rs[0,2], 'ko', label='Initial Pos')

        maxVal = np.max(np.abs(self.rs))

        # ax.set_xlim([-maxVal, maxVal])
        # ax.set_ylim([-maxVal, maxVal])
        # ax.set_zlim([-maxVal, maxVal])

        ax.set_xlabel(['X (km)'])
        ax.set_ylabel(['Y (km)'])
        ax.set_zlabel(['Z (km)'])

        ax.set_aspect('equal')

        ax.set_title(title)

        plt.show()


