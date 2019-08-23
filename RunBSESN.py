from cosmic.evolve import Evolve
import numpy as np
import pandas as pd
import os
call = os.system
import sys
import subprocess
import h5py
import argparse
sys.path.append("./inits")
import run_config
import astropy.units as u
import astropy.constants as c
import BatchRunSN as SN
sys.path.append(run_config.LMXB_INSTALL)

def parse_commandline():
    parser = argparse.ArgumentParser()

    parser.add_argument('--i', type=str, required = True)
    parser.add_argument('--ncores', type = int, default = 16)
    parser.add_argument('--onlyrun', type=str, choices=['bse', 'sn'])
    parser.add_argument('--sn_input', type=str)
    
    return parser.parse_args()

def initialize_binary_table(BSE_input):
    call('python '+run_config.LMXB_INSTALL+'/NewSampling.py --Nsys ' + str(run_config.Nbse) + \
         ' --binaryfile ' + BSE_input + \
         ' --Mexp ' + str(run_config.Mexp))
    
    data = np.loadtxt(BSE_input,skiprows=1).T
    m1, m2, tb, ecc, z, time_to_evolve = data  
    
    tphysf = 10000*np.ones_like(time_to_evolve)
    kstar_1, kstar_2 = np.ones_like(m1),np.ones_like(m2)
    kstar_2[np.where(m2<0.7)[0]] = 0
    columns = ['kstar_1','kstar_2','mass1_binary','mass2_binary','porb','ecc','metallicity','tphysf']
    init_dict = dict(zip(columns,[kstar_1,kstar_2,m1,m2,tb,ecc,z,tphysf]))
    
    return pd.DataFrame(init_dict)
def construct_bse_output(bpp,initC):
    bin_num, evol_type, tphys, Mpre_temp, MdonSN_temp, kstar_Mpre_temp, kstar_MdonSN_temp,epre,Apre = \
                                                                bpp.bin_num.values,bpp.evol_type.values, bpp.tphys.values,\
                                                                bpp.mass_1.values,bpp.mass_2.values, \
                                                                bpp.kstar_1.values,bpp.kstar_2.values,\
                                                                bpp.ecc.values,bpp.sep.values
    
    M1Zams_temp,M2Zams_temp,PZams,eZams =  initC.mass1_binary.values, initC.mass2_binary.values, \
                                       initC.porb.values, initC.ecc.values
    
    i_switch = np.where(evol_type>=16)[0]
    Mpre, MdonSN, kstar_Mpre, kstar_MdonSN = Mpre_temp.copy(), MdonSN_temp.copy(),\
                                             kstar_Mpre_temp.copy(), kstar_MdonSN_temp.copy()
    M1Zams, M2Zams = M1Zams_temp.copy(), M2Zams_temp.copy()
    AZams = period_to_sep(M1Zams,M2Zams,PZams)

    Mpre[i_switch] = MdonSN_temp[i_switch]
    MdonSN[i_switch] = Mpre_temp[i_switch]
    kstar_Mpre[i_switch] = kstar_MdonSN_temp[i_switch]
    kstar_MdonSN[i_switch] = kstar_Mpre_temp[i_switch]
    M1Zams[i_switch] = M2Zams_temp[i_switch]
    M2Zams[i_switch] = M1Zams_temp[i_switch]
    
    cols = ['bin_num','M1Zams','M2Zams','AZams','eZams',\
            'evol_type_preSN','tphys','Mpre','MdonSN','kstar_Mpre','kstar_MdonSN','epre','Apre']
    data = [bin_num, M1Zams, M2Zams, AZams, eZams, \
            evol_type, tphys, Mpre, MdonSN, kstar_Mpre, kstar_MdonSN,epre,Apre]
    return pd.DataFrame(dict(zip(cols,data)))
def get_commit_hashes(run_config):
    batcmd1='git --git-dir='+run_config.COSMIC_INSTALL+'/.git log -1 --format="%H"'
    batcmd2='git --git-dir='+run_config.LMXB_INSTALL+'/.git log -1 --format="%H"'
    
    commit_cosmic = subprocess.check_output(batcmd1, shell=True).strip()
    commit_lmxb = subprocess.check_output(batcmd2, shell=True).strip()
    
    return commit_cosmic, commit_lmxb
def period_to_sep(M1,M2,P):
    G = c.G.value
    M1 = M1*c.M_sun.value
    M2 = M2*c.M_sun.value
    P = P*u.day.to(u.second)
    A = (G*(M1+M2)*4*(P/np.pi)**2)**(1./3.)
    return A*u.m.to(u.Rsun)
if __name__=='__main__':
    args=parse_commandline()
    commit_cosmic, commit_lmxb = get_commit_hashes(run_config)
    
    BSE_input = run_config.BSE_input+args.i+'.in'
    BSE_output = run_config.BSE_output+args.i+'.h5'
    INITC_input = run_config.INITC_input+args.i+'.h5'
    SN_output = run_config.SN_output+args.i+'.h5'

    BSEDict = run_config.BSEDict

    with h5py.File(BSE_output,'w') as f:
        B = f.create_group('BSEDict')

        for key in BSEDict.keys():
            f['BSEDict'][key]=BSEDict[key]
            
        f.attrs['COSMIC_git_hash'] = commit_cosmic
        f.attrs['LMXB_git_hash'] = commit_lmxb

        f.attrs['Nbse'] = run_config.Nbse
        f.attrs['Mexp'] = run_config.Mexp
    
    if args.onlyrun is not None or args.onlyrun == 'bse':

        InitialTable = initialize_binary_table(BSE_input)
        call('rm '+BSE_input)
    
        bpp, bcm, initC  = Evolve.evolve(initialbinarytable=InitialTable, BSEDict=BSEDict,nproc=args.ncores,idx=int(args.i)*run_config.Nbse)
    
        iloc = np.where(((((bpp.evol_type==15) | (bpp.evol_type==15.5)) & (bpp.kstar_2<13)) |
            (((bpp.evol_type==16) | (bpp.evol_type==16.5)) & (bpp.kstar_1<13)))&(bpp.sep>0))[0]
        bpp_preSN = bpp.iloc[iloc]
    
        bin_num=bpp.iloc[iloc].index
        initC = initC.set_index('bin_num')
        initC_preSN = initC.loc[bin_num]
        initC_preSN.to_hdf(INITC_input,key='initC')
    
        preSNbinaries = construct_bse_output(bpp_preSN,initC_preSN)
    
        preSNbinaries.to_hdf(BSE_output,key='bse')
    
    if args.onlyrun is not None:
        call("sed '9i#SBATCH --error=output/" + i + "arrayJob_%A_%a.err' SNsubmit.sh > " + i + "_SNsubmit.sh")
        call("sed -i '10i#SBATCH --output=output/" + i + "arrayJob_%A_%a.out' " + i + "_SNsubmit.sh")
        with open(i + "_SNsubmit.sh", a) as f:
            f.write('python BatchRunSN.py -x ${SLURM_ARRAY_TASK_ID} --nkicks ' + run_config.Nkick + ' --BSEnum' + i + '--rows ' + 1000)
        #call('sbatch ' + i + '_SNsubmit.sh')
        
    elif args.onlyrun == 'sn':
        if args.sn_input is None:
            parser.error('sn_input is required when onlyrun = sn')
        else:
            postSN = SN.runSN(args.i, args.sn_input, 100, 1000)
            postSN.to_hdf(SN_output,key='sn')
