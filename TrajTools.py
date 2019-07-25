from astropy import units as u
import numpy as np
gauss=np.random.normal
import scipy.stats as stats
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy import constants as C
Msun=C.M_sun.value
m_kpc=u.kpc.to(u.m)

G=C.G.value
s_year=31556952.0

def potential(x,y,z):
    x  = x*u.kpc.to(u.m)
    y  = y*u.kpc.to(u.m)
    z  = z*u.kpc.to(u.m)
    Md = 1.45e11*Msun
    Mb = 9.3e9*Msun
    Mn = 1e10*Msun
    bd = 5.5*u.kpc.to(u.m)
    bb = 0.25*u.kpc.to(u.m)
    bn = 1.5*u.kpc.to(u.m)
    h1 = 0.325*u.kpc.to(u.m)
    h2 = 0.090*u.kpc.to(u.m)
    h3 = 0.125*u.kpc.to(u.m)
    aG = 2.4*u.kpc.to(u.m)
    B1 = 0.4
    B2 = 0.5
    B3 = 0.1
    S1 = B1*(z**2+h1**2)**0.5
    S2 = B2*(z**2+h2**2)**0.5
    S3 = B3*(z**2+h3**2)**0.5
    vd = (-G*Md)/((aG+S1+S2+S3)**2+bd**2+x**2+y**2)**0.5
    vb = (-G*Mb)/(bb**2+x**2+y**2)**0.5
    vn = (-G*Mn)/(bn**2+x**2+y**2)**0.5
    vg = vn+vb+vd
    return vg

def vdot(t,R):
        M=[(1.45e11)*Msun,(9.3e9)*Msun,(1.0e10)*Msun]
        B=[0.4,0.5,0.1]
        H=[0.325*m_kpc,0.090*m_kpc,0.125*m_kpc]
        Ag=2.4*m_kpc
        b=[5.5*m_kpc,0.25*m_kpc,1.5*m_kpc]
        x,y,z=R[3],R[4],R[5]
        Tx0=-(G*M[0]*x*pow((b[0]**2)+(x**2)+(y**2)+((Ag+(B[0]*pow((H[0]**2)+(z**2),0.5))+(B[1]*pow((H[1]**2)+(z**2),0.5))+(B[2]*pow((H[2]**2)+(z**2),0.5)))**2),-3.0/2))
        Tx1=-(G*M[1]*x*pow((b[1]**2)+(x**2)+(y**2),-3.0/2))
        Tx2=-(G*M[2]*x*pow((b[2]**2)+(x**2)+(y**2),-3.0/2))

        Ty0=-(G*M[0]*y*pow((b[0]**2)+(x**2)+(y**2)+((Ag+(B[0]*pow((H[0]**2)+(z**2),0.5))+(B[1]*pow((H[1]**2)+(z**2),0.5))+(B[2]*pow((H[2]**2)+(z**2),0.5)))**2),-3.0/2))
        Ty1=-(G*M[1]*y*pow((b[1]**2)+(x**2)+(y**2),-3.0/2))
        Ty2=-(G*M[2]*y*pow((b[2]**2)+(x**2)+(y**2),-3.0/2))

        h0term=pow((H[0]**2)+(z**2),0.5)
        h1term=pow((H[1]**2)+(z**2),0.5)
        h2term=pow((H[2]**2)+(z**2),0.5)
        Tz1=-((B[0]*z/h0term)+(B[1]*z/h1term)+(B[2]*z/h2term))
        Tz2=(Ag+(B[0]*h0term)+(B[1]*h1term)+(B[2]*h2term))
        Tz3=pow((b[0]**2)+(x**2)+(y**2)+(Tz2**2),-3.0/2)
        return np.array([(Tx0+Tx1+Tx2),(Ty0+Ty1+Ty2),(G*M[0]*Tz1*Tz2*Tz3),R[0],R[1],R[2]])

def Vrad(X01,Y01,U1,V1):
      r=np.sqrt((X01**2)+(Y01**2))
      return ((U1*X01)+(V1*Y01))/r
def Vcirc(U,V,W,Vr):
      v2=(V**2)+(U**2)
      
      return np.sqrt(v2-(Vr**2))
def getPec(X01,Y01,Z01,U,V,W):
      vrad=Vrad(X01,Y01,U,V)
      vcirc=Vcirc(U,V,W,vrad)
      return np.sqrt((vrad**2)+((vcirc-getVrot(X01,Y01,Z01))**2)+(W**2))
    
def getPecFixed(X01,Y01,Z01,U,V,W):
      vrad=Vrad(X01,Y01,U,V)      
      vcirc=Vcirc(U,V,W,vrad)
      return np.sqrt((vrad**2)+((vcirc-238000.)**2)+(W**2))

def getVrot(X,Y,Z): #m
    R=[0,0,0,X,Y,0]
    P0=-1*((X*vdot(0,R)[0])+(Y*vdot(0,R)[1]))
    Vrot=np.sqrt(P0)#*numpy.sqrt((X**2)+(Y**2)))
    return Vrot   #m/s

Rs = 8.05 #kpc Miller Jones
Omega = getVrot(-Rs*u.kpc.to(u.m),0,0)*u.m.to(u.km)
pmsun=[11.1,12.24+Omega,7.25] #km/s Miller Jones
print ('Omega = '+str(Omega))
def getnonPec(X,Y,Z,Up,Vp,Wp):
    Vrot=getVrot(X,Y,Z)
    R=np.sqrt((X**2)+(Y**2))
    U=Up+(Vrot*Y/R)
    V=Up+(Vrot*(-X/R))
    return U,V,Wp
def drawGauss(ARGS):

    ARGS_RAND=[gauss(arg[0],arg[1]) for arg in ARGS]
    return ARGS_RAND
def getRandXYZUVW(ra,dec,distance,pm_ra,pm_dec,radial_velocity,v_sun=pmsun,galcen_distance=Rs,dlow=None,dhigh=None,d_musig=None):
        PM_RA,PM_DEC,RADIAL_VELOCITY=drawGauss([pm_ra,pm_dec,radial_velocity])
        if d_musig:
            DISTANCE=0.0
            while DISTANCE<dlow or DISTANCE>dhigh:
                DISTANCE = stats.lognorm(s = d_musig[1],scale = np.exp(d_musig[0])).rvs()
        else:
            DISTANCE = drawGauss([distance])[0]
        
        return getXYZUVW(ra,dec,DISTANCE,PM_RA,PM_DEC,RADIAL_VELOCITY,v_sun=v_sun,galcen_distance=galcen_distance) #kpc and km/s
def getXYZUVW(ra,dec,distance,pm_ra_cosdec,pm_dec,radial_velocity,v_sun=pmsun,galcen_distance=Rs):
    #degree,degree,kpc,mas/yr,mas/yr,km/s,km/s,kpc
    c1=coord.ICRS(ra=ra*u.degree,dec=dec*u.degree,distance=distance*u.kpc,\
        pm_ra_cosdec=pm_ra_cosdec*u.mas/u.yr,pm_dec=pm_dec*u.mas/u.yr,\
        radial_velocity=radial_velocity*u.km/u.s)

    gc_frame=coord.Galactocentric(galcen_distance=galcen_distance*u.kpc,\
        galcen_v_sun=coord.CartesianDifferential(v_sun*u.km/u.s))
    gc2=c1.transform_to(gc_frame)
    #kpc,kpc,kpc,km/s,km/s,km/s
    return [gc2.x.value,gc2.y.value,gc2.z.value,gc2.v_x.value,gc2.v_y.value,gc2.v_z.value]








#Xsys,Ysys,Zsys,Usys,Vsys,Wsys=getXYZUVW(ra[0],dec[0],distance[0],pm_ra[0],pm_dec[0],radial_velocity[0])
#print Usys,Vsys,Wsys



