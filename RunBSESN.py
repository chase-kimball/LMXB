from cosmic.evolve import Evolve
import numpy as np
import pandas as pd
import os
call = os.system
import sys
import subprocess
import h5py
sys.path.append("./inits")
import run_config
sys.path.append(run_config.LMXB_INSTALL)

def parse_commandline():
    
    parser.add_argument('--i', type=str, required = True)
    parser.add_argument('--ncores', type = int, default = 16)

    return parser.parse_args()

def initialize_binary_table(BSE_input):
    call('python NewSampling.py --Nsys ' + str(run_config.Nbse) + \
         ' --binaryfile ' + run_config.BSE_input + args.i + \
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
    evol_type, tphys, Mpre_temp, MdonSN_temp, kstar_Mpre_temp, kstar_MdonSN_temp,epre,Apre = \
                                                                bpp.evol_type.values, bpp.tphys.values,\
                                                                bpp.mass_1.values,bpp.mass_2.values, \
                                                                bpp.kstar_1.values,bpp.kstar_2.values,\
                                                                bpp.ecc.values,bpp.sep.values
    
    index, M1Zams_temp,M2Zams_temp,PZams,eZams = initC.bin_num.values, initC.mass_1.values, initC.mass_2.values \
                                       initC.porb.values, initC.ecc.values
    AZams = period_to_sep(PZams)
    
    i_switch = np.where(evol_type>=16)[0]
    Mpre, MdonSN, kstar_Mpre, kstar_MdonSN = Mpre_temp.copy(), MdonSN_temp.copy(),\
                                             kstar_Mpre_temp.copy(), kstar_MdonSN_temp.copy()
    M1Zams, M2Zams = M1Zams_temp.copy(), M2Zams_temp.copy()
    
    Mpre[i_switch] = MdonSN_temp[i_switch]
    MdonSN[i_switch] = Mpre_temp[i_switch]
    kstar_Mpre[i_switch] = kstar_MdonSN_temp[i_switch]
    kstar_MdonSN[i_switch] = kstar_Mpre_temp[i_switch]
    M1Zams[i_switch] = M2Zams_temp[i_switch]
    M2Zams[i_switch] = M1Zams_temp[i_switch]
    
    cols = ['bin_num','M1Zams','M2Zams','AZams','eZams',\
            'evol_type_preSN','tphys','Mpre','MdonSN','kstar_Mpre','kstar_MdonSN','epre','Apre']
    data = [index, M1Zams, M2Zams, AZams, eZams, \
            evol_type, tphys, Mpre, MdonSN, kstar_Mpre, kstar_MdonSN,epre,Apre]
    return pd.DataFrame(dict(zip(cols,data)))
def get_commit_hashes(run_config):
    batcmd1='git --git-dir='+run_config.COSMIC_INSTALL+'/.git log -1 --format="%H"'
    batcmd2='git --git-dir='+run_config.LMXB_INSTALL+'/.git log -1 --format="%H"'
    
    commit_cosmic = subprocess.check_output(batcmd1, shell=True).strip()
    commit_lmxb = subprocess.check_output(batcmd2, shell=True).strip()
    
    return commit_cosmic, commit_lmxb

if __name__=='__main__':
    args=parse_commandline()
    commit_cosmic, commit_lmxb = get_commit_hashes(run_config)
    
    BSE_input = run_config.BSE_input+args.i+'.in'
    BSE_output = runconfig.BSE_output+args.i+'.h5'
    INITC_input = run_config.INITC_input+args.i+'.h5'
    
    BSEdict = run_config.BSEdict

    with h5py.File(BSE_ouput,'w') as f:
        B = f.create_group('BSEdict')

        for key in BSEDict.keys():
            f['BSEdict'][key]=BSEDict[key]
            
        f.attrs['COSMIC_git_hash'] = commit_cosmic
        f.attrs['LMXB_git_hash'] = commit_lmxb

        f.attrs['Nbse'] = run_config.Nbse
        f.attrs['Mexp'] = run_config.Mexp
    


    
    
    InitialTable = initialize_binary_table(BSE_input)
    call('rm '+BSE_input)
    

    bpp, bcm, initC  = Evolve.evolve(initialbinarytable=InitTable, BSEDict=BSEDict,nproc=args.ncores,args.i*run_config.Nbse)
    
    iloc = np.where(((((bpp.evol_type==15) | (bpp.evol_type==15.5)) & (bpp.kstar_2<13)) |
        (((bpp.evol_type==16) | (bpp.evol_type==16.5)) & (bpp.kstar_1<13)))&(bpp.sep>0))[0]
    bpp_preSN = bpp[iloc]
    
    bin_num=bpp.iloc[iloc].index
    initC = initC.set_index('bin_num')
    initC_preSN = initC.loc[bin_num]
    initC_preSN.to_hdf(INITC_input,key='initC')
    
    preSNbinaries = construct_bse_output(bpp_preSN,initC_preSN)
    
    
    preSNbinaries.to_hdf(BSE_output,key='bse')
    
