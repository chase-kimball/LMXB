import TrajTools
import astropy.units as u
def getSystem(sysname):
    Sys_Dict ={

        'GRS1915': {
            'DX' : [0.5, 0.2, 0.5, 0.025],
            'Mbh': [12.4,2.0,1.8],
            'Mbh_musig':[2.5058656370,0.153823506514],
            'q': [0.042,0.024],
            'P' : [33.85, 0.0],
            'Teff':[4766.5,333.25],
            'ra': 288.7981250,
            'dec' : 10.9457667,
            'd' : [8.6,2.0, 1.6],
            'pm_ra' : [-3.19,0.03],#actually pm_ra*cos(dec)
            'pm_dec':[-6.24,0.05],
            'v_rad':[12.3,1.0],
            'dlow':5.9,
            'dhigh':12.0,
            'd_musig':[2.13026759674,0.207338402227]


        },
        
        'V404': {
            'DX': [0.5, 0.1, 0.5, 0.025],
            'Mbh':[11.5,1.75],
            'q_musig':[-2.81627856718, 0.0757344098565],
            'q': [0.06,0.004,.005],
            'P': [6.4714, 0.0],
            'Teff':[4400, 100.0],
            'ra': 306.0158333,
            'dec': 33.8675556,
            'd': [2.39,0.14],
            'pm_ra': [-5.04,0.22],
            'pm_dec': [-7.64,0.03],
            'v_rad': [-0.4,2.2],
            'figsize':(20,11),

            'binZAMS':[20,15,50,15],
            'bwZAMS':[2,.05,500,.05],

            'binPre':[15,15,50,50],
            'bwPre':[.25,.1,1.5,15],

            'binPost':[15,10,90,100],
            'bwPost':[.1,.1,5.0,.05],

            'binRLO':[[4.6,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5],12,[0.5,1.3,2.1,2.9]],
            'bwRLO':[0.5,0.5,0.4],

            'fullname': 'V404 Cygni',
            'dlow':None,
            'dhigh':None,
            'd_musig':None
        }
    }

    system = Sys_Dict[sysname]
    Xsys,Ysys,Zsys,Usys,Vsys,Wsys = TrajTools.getXYZUVW(system['ra'],system['dec'],system['d'][0],system['pm_ra'][0],system['pm_dec'][0],system['v_rad'][0])
    system['Xsys'],system['Ysys'],system['Zsys'],system['Usys'],system['Vsys'],system['Wsys']=Xsys,Ysys,Zsys,Usys,Vsys,Wsys
    system['XYZUVW'] = [Xsys,Ysys,Zsys,Usys,Vsys,Wsys]

    print# (TrajTools.getPec(Xsys*u.kpc.to(u.m),Ysys*u.kpc.to(u.m),Zsys*u.kpc.to(u.m),Usys*1e3,Vsys*1e3,Wsys*1e3)*1e-3)
    return system
    
        
