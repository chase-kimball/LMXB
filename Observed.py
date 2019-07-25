def getSysDict(sysname):
    Sys_Dict ={

        'GRS1915': {
            'ra': 288.7981250,
            'dec' : 10.9457667,
            'd' : [8.6,2.0, 1.6],
            'pm_ra' : [-3.19,0.03],#actually pm_ra*cos(dec)
            'pm_dec':[-6.24,0.05],
            'v_rad':[12.3,1.0],
            
            'fullname': 'GRS 1915+105'

        },
        
        'V404': {
            'ra': 306.0158333,
            'dec': 33.8675556,
            'd': [2.39,0.14],
            'pm_ra': [-5.04,0.22],
            'pm_dec': [-7.64,0.03],
            'v_rad': [-0.4,2.2],

            'fullname': 'V404 Cygni'
        },
        
        'SCO': {
            'l': 359.09,
            'b': 23.78,
            'd': [2.8,0.3],
            'pm_ra': [-6.88,0.07],
            'pm_dec': [-12.02,0.16],
            'v_rad': [-113.8,0.6],
            
            'fullname': 'SCO X-1'
        },
        
        'LS': {
            'l': 16.88,
            'b': -1.29,
            'd': [2.9,0.3],
            'pm_ra': [7.10,0.13],
            'pm_dec': [-8.75,0.16],
            'v_rad': [17.2,0.5],
            
            'fullname': 'LS 5039'
        },
        
        'Aql': {
            'l': 34.67,
            'b': -4.68,
            'd': [5.0,0.9],
            'pm_ra': [-2.64,0.14],
            'pm_dec': [-3.53,1.40],
            'v_rad': [30,10],
            
            'fullname': 'Aql X-1'
        },
        
        'Cyg': {
            'l': 87.33,
            'b': -11.32,
            'd': [11,2],
            'pm_ra': [-3.00,0.68],
            'pm_dec': [-0.64,0.68],
            'v_rad': [-209.6,0.8],
            
            'fullname': 'Cyg X-2'
        },
        
        'XTEJ1118':{
            'l': 157.66,
            'b': 62.32,
            'd': [1.72, 0.10],
            'pm_ra': [-16.8, 1.6],
            'pm_dec': [-7.4, 1.6],
            'v_rad': [2.7, 1.1],
            
            'fullname': 'XTE J1118+480'
        },
        
        'GROJ1655':{
            'l': 344.98,
            'b': 2.46,
            'd': [3.2, 0.2],
            'pm_ra': [-3.3, 0.50],
            'pm_dec': [-4.0, 0.4],
            'v_rad': [-141.9, 1.3],
            
            'fullname': 'GRO J1655-40'
        },
        
        'LSI61':{
            'l': 135.68,
            'b': 1.09,
            'd': [2.0, 0.2],
            'pm_ra': [-0.30, 0.07],
            'pm_dec': [-0.26, 0.05],
            'v_rad': [-40.2, 1.9],
            
            'fullname': 'LSI +61 303'
        },
        
        'CEN':{
            'l': 332.24,
            'b': 23.88,
            'd': [1.4, 0.3],
            'pm_ra': [-56, 10],
            'pm_dec': [11, 10],
            'v_rad': [189.6, 0.2],
            
            'fullname': 'Cen X-4'
        }
    }

    system = Sys_Dict[sysname]
    return system

import Observed #dictionary file
import TrajTools
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import pdb

gauss=np.random.normal

import scipy.stats as stats
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy import constants as C

Msun=C.M_sun.value
m_kpc=u.kpc.to(u.m)


## Class for systems

class LMXB_Sys:
    def __init__(self, sysname):
        self.system = Observed.getSysDict(sysname)
        
        ## Assigns the parameters of each system to a variable
        self.l = self.system.get('l')
        self.b = self.system.get('b')
        self.d = self.system.get('d')
        self.pm_ra = self.system.get('pm_ra')
        self.pm_dec = self.system.get('pm_dec')
        self.v_rad = self.system.get('v_rad')
        self.fullname = self.system.get('fullname')
        
        ## Converts into earth centric coordinates
        Sys_gal = SkyCoord(l = self.system.get('l')*u.deg, b = self.system.get('b')*u.deg, frame = 'galactic')
        self.ra = Sys_gal.icrs.ra.value
        self.dec = Sys_gal.icrs.dec.value
             
        def getCartesian(self):
            self.X0, self.Y0, self.Z0, self.U0, self.V0, self.W0 = TrajTools.getXYZUVW(self.ra, self.dec, self.d[0], self.pm_ra[0], self.pm_dec[0], self.v_rad[0])
        getCartesian(self)
        
    def setRandUVWXYZ(self):
        self.X0, self.Y0, self.Z0, self.U0, self.V0, self.W0 = TrajTools.getRandXYZUVW(self.ra, self.dec, self.d, self.pm_ra, self.pm_dec, self.v_rad)
        
        
    # Numerical Integrator    
    def euler(self, R0, dt, T, force_function):
        time = 0
        tarr = [time]
        R = R0
        print(R)
        #pdb.set_trace()
        R_t = [R0] 
        while np.abs(time) < T:
            R = R + dt*force_function(time, R)
            time += dt
            R_t.append(R)
            tarr.append(time)
        R_t = np.array(R_t)

        return R_t.T, np.array(tarr)
    
    def RK4(self, R0, dt, T, force_function):
        t = 0
        tarr = [t]
        R = R0
        R_t = [R0] 
        while np.abs(t) < T:
            R1 = R
            t1 = t
            vdot1 = force_function(t1, R1)

            R2 = R + 0.5*dt*vdot1
            t2 = t + 0.5*dt
            vdot2 = force_function(t2, R2)


            R3 = R + 0.5*dt*vdot2
            t3 = t + 0.5*dt
            vdot3 = force_function(t3, R3)


            R4 = R + dt*vdot3
            t4 = t + dt
            vdot4 = force_function(t4, R4)

            R = R + (dt/6)*(vdot1 + 2*vdot2 + 2*vdot3 + vdot4)
            t = t + dt

            R_t.append(R)
            tarr.append(t)

        R_t = np.array(R_t)

        return R_t.T, np.array(tarr)

    
        #def gravity(R, t):
            #g = 9.8
            #return np.array([-g, R[0]])
        
        
    def plotTrajectory(self, T, dt, integrator):
        
        R0 = np.array([self.U0*1000, self.V0*1000, self.W0*1000, self.X0*u.kpc.to(u.m), self.Y0*u.kpc.to(u.m), self.Z0*u.kpc.to(u.m)]) 
        #print(R0,dt*u.Gyr.to(u.s),T*u.Gyr.to(u.s), TrajTools.vdot(T, R0))
        R_t, t = integrator(R0,dt*u.Gyr.to(u.s),T*u.Gyr.to(u.s), TrajTools.vdot)
        self.R_t = R_t
        self.R_t[0] = R_t[0]/1000
        self.R_t[1] = R_t[1]/1000
        self.R_t[2] = R_t[2]/1000
        self.R_t[3] = R_t[3]*u.m.to(u.kpc)
        self.R_t[4] = R_t[4]*u.m.to(u.kpc)
        self.R_t[5] = R_t[5]*u.m.to(u.kpc)
        self.t = t
        print(self.U0,self.V0,self.W0,self.X0,self.Y0,self.Z0)
        print(R0)
        plt.plot(self.R_t[3], self.R_t[4])
    
       
    def getTrajectory(self, T, dt, integrator):
        
        R0 = np.array([self.U0*1000, self.V0*1000, self.W0*1000, self.X0*u.kpc.to(u.m), self.Y0*u.kpc.to(u.m), self.Z0*u.kpc.to(u.m)]) 
        #print(R0,dt*u.Gyr.to(u.s),T*u.Gyr.to(u.s), TrajTools.vdot(T, R0))
        R_t, t = integrator(R0,dt*u.Gyr.to(u.s),T*u.Gyr.to(u.s), TrajTools.vdot)
        self.R_t = R_t
        self.R_t[0] = R_t[0]/1000
        self.R_t[1] = R_t[1]/1000
        self.R_t[2] = R_t[2]/1000
        self.R_t[3] = R_t[3]*u.m.to(u.kpc)
        self.R_t[4] = R_t[4]*u.m.to(u.kpc)
        self.R_t[5] = R_t[5]*u.m.to(u.kpc)
        self.t = t
        print(R_t)
        return R_t, t 
        
        
                