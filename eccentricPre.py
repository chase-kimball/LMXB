# -*- coding: utf-8 -*-
# Copyright (C) Charles Kimall and Michael Zevin (2018)
#
# This file is part of the progenitor package.
#
# progenitor is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# progenitor is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with progenitor.  If not, see <http://www.gnu.org/licenses/>.

__author__ = ['Michael Zevin <michael.zevin@ligo.org>', 'Chase Kimball <charles.kimball@ligo.org']
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['System']


import numpy as np
import pandas as pd

import astropy.units as units
import astropy.constants as constants

from scipy.integrate import ode
from scipy.integrate import quad
from scipy.optimize import brentq
import zams as zams

class System:
    """
    Places system described by Mhe, Mcomp, Apre, epre and position r(R,galphi,galcosth) in galaxy model gal

    Applies SNkick Vkick and mass loss Mhe-Mns to obtain Apost, epost, and SN-imparted systemic velocity V
    
    """
    def __init__(self, Mcomp, Mhe, Apre, epre, Nkick=1000, Vkick=None, Mns=None, sys_flag=None, galphi=None, galcosth=None, omega=None, phi=None, costh=None,th_ma = None):
        """ 
        #Masses in Msun, Apre in Rsun, Vkick in km/s, R in kpc
        #galphi,galcosth,omega, phi, costh (position, initial velocity, and kick angles) sampled randomly, unless specified (>-1)
        #galphi, galcosth correspond to azimuthal and polar angles -- respectively --  in the galactic frame
        #phi, costh are defined in comments of SN:
        #   theta: angle between preSN He core velocity relative to Mcomp (i.e. the positive y axis) and the kick velocity
        #   phi: angle between Z axis and projection of kick onto X-Z plane
        #omega: angle between the galactic velocity corresponding to a circular orbit in the r-z plane and
        #the actual galactic velocity preSN corresponding to a circular orbit
        
        """
    
        # Convert inputs to SI


        self.sys_flag = sys_flag
        self.Nkick = Nkick

        if np.any(Vkick): self.Vkick = Vkick*units.km.to(units.m)
        else: self.Vkick = np.random.uniform(0,1000,self.Nkick)*units.km.to(units.m)
        if np.any(phi): self.phi = phi
        else: self.phi = np.random.uniform(0,2*np.pi,self.Nkick)

        if np.any(costh): self.costh = costh
        else: self.costh = np.random.uniform(-1,1,self.Nkick)
        if np.any(Mns): self.Mns = Mns*units.M_sun.to(units.kg)
        else: self.Mns = np.random.uniform(3.,Mhe,self.Nkick)*units.M_sun.to(units.kg)
            
        if np.any(th_ma): self.th_ma = th_ma
        else: self.th_ma = np.random.uniform(0,2*np.pi,self.Nkick)
        self.E_ma =np.array([brentq(lambda x:ma -x + epre*np.sin(x),0,2*np.pi) for ma in self.th_ma])
        self.rpre = Apre*(1.-epre*np.cos(self.E_ma))*units.R_sun.to(units.m)
        self.Mhe = np.full((self.Nkick,), Mhe)*units.M_sun.to(units.kg)
        self.Mcomp = np.full((self.Nkick,), Mcomp)*units.M_sun.to(units.kg)
        self.Apre = np.full((self.Nkick,),Apre)*units.R_sun.to(units.m)
        self.epre = np.full((self.Nkick,),epre)
        
        # Get projection of R in the x-y plane to save later into output file

    def SN(self):
        """
        
        Mhe lies on origin moving in direction of positive y axis, Mcomp on negative X axis, Z completes right-handed coordinate system
        
        theta: angle between preSN He core velocity relative to Mcomp (i.e. the positive y axis) and the kick velocity
        phi: angle between Z axis and projection of kick onto X-Z plane
        
        Vr is velocity of preSN He core relative to Mcomp, directed along the positive y axis

        Vkick is kick velocity with components Vkx, Vky, Vkz in the above coordinate system
        V_sys is the resulting center of mass velocity of the system IN THE TRANSLATED COM FRAME, imparted by the SN

        Paper reference:

        Kalogera 1996: http://iopscience.iop.org/article/10.1086/177974/meta
            We use Eq 1, 3, 4, and 34: giving Vr, Apost, epost, and (Vsx,Vsy,Vsz) respectively
            Also see Fig 1 in that paper for coordinate system

        
        """

        self.flag=0      # set standard flag        

        G = constants.G.value
        Mhe, Mcomp, Mns, Apre, epre, rpre, Vkick, costh, phi = self.Mhe, self.Mcomp, self.Mns, self.Apre, self.epre, self.rpre, self.Vkick, self.costh, self.phi


        sinth = np.sqrt(1-(costh**2))
        #Mhe lies on origin moving in direction of positive y axis, Mcomp on negative X axis, Z completes right-handed coordinate system
        #See Fig 1 in Kalogera 1996

        # theta: angle between preSN He core velocity relative to Mcomp (i.e. the positive y axis) and the kick velocity
        # phi: angle between Z axis and projection of kick onto X-Z plane
        Vkx = Vkick*sinth*np.sin(phi)
        Vky = Vkick*costh
        Vkz = Vkick*sinth*np.cos(phi)
        
        #Eq 1, Kalogera 1996
        Vr = np.sqrt(G*(Mhe+Mcomp)*(2./rpre - 1./Apre))
        Mtot=Mns+Mcomp

        #Eqs 3 and 4, Kalogera 1996
        Apost = ((2.0/rpre) - (((Vkick**2)+(Vr**2)+(2*Vky*Vr))/(G*Mtot)))**-1
        sin_ang_ma = np.round(np.sqrt(G*(Mhe+Mcomp)*(1-epre**2)*Apre)/(rpre*Vr),5)
        self.sin_ang_ma = sin_ang_ma
        cos_ang_ma = np.sqrt(1-sin_ang_ma**2)
        self.cos_ang_ma = cos_ang_ma
        x = ((Vkz**2)+((sin_ang_ma*(Vr+Vky)-cos_ang_ma*(Vkx))**2))*(rpre**2)/(G*Mtot*Apost)
        epost = np.sqrt(1-x)
        # Eq 34, Kalogera 1996
        VSx = Mns*Vkx/Mtot
        VSy = (1.0/Mtot)*((Mns*Vky)-((Mhe-Mns)*Mcomp*Vr/(Mhe+Mcomp)))
        VSz = Mns*Vkz/Mtot
        V_sys = np.sqrt((VSx**2)+(VSy**2)+(VSz**2))

        self.Apost, self.epost, self.VSx, self.VSy, self.VSz, self.V_sys, self.Vr = Apost, epost, VSx, VSy, VSz, V_sys, Vr
        
        def SNCheck(self):
            """
            Paper References:

            Willems et al 2002: http://iopscience.iop.org/article/10.1086/429557/meta
                We use eq 21, 22, 23, 24, 25, 26 for checks of SN survival

            Kalogera and Lorimer 2000: http://iopscience.iop.org/article/10.1086/308417/meta

            
            V_He;preSN is the same variable as V_r from Kalogera 1996
            
            """
            Mhe, Mcomp, Mns, Apre, epre, rpre, Apost, epost, Vr, Vkick = self.Mhe, self.Mcomp, self.Mns, self.Apre, self.epre, self.rpre, self.Apost, self.epost, self.Vr, self.Vkick
            #Equation numbers and quotes in comments correspond to Willems et al. 2002 paper on J1655.
            Mtot_pre = Mhe + Mcomp
            Mtot_post = Mns + Mcomp

            # SNflag1: eq 21 (with typo fixed). Continuity demands Post SN orbit must pass through preSN positions.
            #from Flannery & Van Heuvel 1975                                                             

            self.oldSNflag1 = (1-epost <= Apre*(1-epre)/Apost) & (Apre*(1+epre)/Apost <= 1+epost)
            self.SNflag1 = (1-epost <= rpre/Apost) & (rpre/Apost <= 1+epost)

            # SNflag2: Equations 22 & 23. "Lower and upper limits on amount of orbital contraction or expansion that can take place                                
            #for a given amount of mass loss and a given magnitude of the kick velocity (see, e.g., Kalogera & Lorimer 2000)"                         self.Mhe, self.Mcomp,self.Apre, self.epre = Mhe, Mcomp, Apre, epre           

            self.SNflag2 = (Apre/Apost < 2-((Mtot_pre/Mtot_post)*((Vkick/Vr)-1)**2)) & (Apre/Apost > 2-((Mtot_pre/Mtot_post)*((Vkick/Vr)+1)**2))

            #SNflag3: Equations 24 and 25."The magnitude of the kick velocity imparted to the BH at birth is restricted to the
            #range determined by (Brandt & Podsiadlowski 1995; Kalogera & Lorimer 2000)
            #the first inequality expresses the requirement that the binary must remain bound after the SN explosion,
            #while the second inequality yields the minimum kick velocity required to keep the system bound if more than
            #half of the total system mass is lost in the explosion.

            self.SNflag3 = (Vkick/Vr < 1 + np.sqrt(2*Mtot_post/Mtot_pre)) & ((Mtot_post/Mtot_pre > 0.5) | (Vkick/Vr>1 - np.sqrt(2*Mtot_post/Mtot_pre)))
            
            

            #SNflag4: Eq 26 "An upper limit on the mass of the BH progenitor can be derived from the condition that the
            #azimuthal direction of the kick is real (Fryer & Kalogera 1997)"
            self.SNflag4 = np.full(self.SNflag3.shape,0)
            ieGT1 = np.where((epost>1)|(epost<0))[0]
            ieLT1 = np.where(epost<=1)[0]
            self.SNflag4[ieGT1] = False
            kvar=2*(Apost[ieLT1]/Apre[ieLT1])-(((Vkick[ieLT1]**2)*Apost[ieLT1]/(G*Mtot_post[ieLT1]))+1)

            tmp1 = kvar**2 * Mtot_post[ieLT1] * (Apre[ieLT1]/Apost[ieLT1])
            tmp2 = 2 * (Apost[ieLT1]/Apre[ieLT1])**2 * (1-epost[ieLT1]**2) - kvar
            tmp3 = - 2 * (Apost[ieLT1]/Apre[ieLT1]) * np.sqrt(1-epost[ieLT1]**2) * np.sqrt((Apost[ieLT1]/Apre[ieLT1])**2 * (1-epost[ieLT1]**2) - kvar)
            prgmax = -Mcomp[ieLT1] + tmp1 / (tmp2 + tmp3)

            self.SNflag4[ieLT1] = Mhe[ieLT1] <= prgmax
            # FIX ME: additionally, Kalogera 1996 mentions requirement that NS stars don't collide
            # Apost*(1-epost)> Rns1+Rns2    (eq 16 in that paper)
            # Is there analytic expression for NS radius?
            self.ieGT1 = ieGT1
            self.ieLT1 = ieLT1
            self.SNflag5 = (1/epost)*(1.-rpre/Apost)<=1.
            
            self.SNflag6 = zams.rzams(Mcomp)<zams.roche_lobe(Mcomp,Mns)*Apost
             
            self.SNflag7 = zams.rhe(Mhe)<zams.roche_lobe(Mhe,Mcomp)*Apost
            
            self.SNflags = [self.SNflag1, self.SNflag2, self.SNflag3, self.SNflag4, self.SNflag5, self.SNflag6, self.SNflag7]

        SNCheck(self)
        #self.Mns, self.Mcomp, self.Mhe = np.asarray([self.Mns,self.Mcomp,self.Mhe])*u.kg.to(u.M_sun)
        #self.Vkick = self.Vkick*u.km.to(u.m)
        #self.Apre = self.Apre*u.m.to(u.R_sun)
        #self.Apost = self.Apost*u.m.to(u.R_sun)

    



