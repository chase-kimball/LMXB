import random
import numpy as np
import os
import argparse
import decimal
import warnings
import sys
sys.path.append('..')
from multiprocessing import Pool
from functools import partial
from operator import is_not
from scipy import optimize
import logging
import astropy.units as u
import SNReviewed as SN
def parse_commandline():
        parser = argparse.ArgumentParser()
        
        parser.add_argument('--SNoutfile', type = str, required = True)
        parser.add_argument('--BSEoutfile', type = str, required = True)
        parser.add_argument('--NSNbatch', type = int, default = 20)
        parser.add_argument('--ncores', type = int, default = 16)
        args = parser.parse_args()

        return args
def deploy(bse):
    Nkick = 1000

    Mpre0=bse[6]                  ##### BH progenitor
    Mdons0=bse[7]                          ##### Donor pre-SN
    Apre0=bse[8]                           ##### Orbital sep pre-SN 
    epre0=bse[9]

    
    SS = SN.System(Mdons0,Mpre0,Apre0,epre0,Nkick=Nkick)
    SS.SN()
    
    iS = np.where(SS.SNflag1*SS.SNflag2*SS.SNflag3*SS.SNflag4*SS.SNflag5*SS.SNflag6*SS.SNflag7==1)[0]
    print iS
    shape = (len(iS),)
    
    index=np.full(shape,bse[0])
    M1Zams=np.full(shape,bse[1])
    M2Zams=np.full(shape,bse[2])
    AZams=np.full(shape,bse[4])
    eZams=np.full(shape,bse[5])
    Mpre=np.full(shape,Mpre0)                  ##### BH progenitor
    Mdons=np.full(shape,Mdons0)                          ##### Donor pre-SN
    Apre=np.full(shape,Apre0)                           ##### Orbital sep pre-SN 
    epre=np.full(shape,epre0)
    thatflag=np.full(shape,bse[-2])
    switchflag=np.full(shape,bse[-6])
    disdum=np.full(shape,bse[-7])
    codum=np.full(shape,bse[-8])
    cedum=np.full(shape,bse[13])
    
    prog = np.array([index,M1Zams,M2Zams,AZams,eZams,Mpre,Mdons,Apre,epre,SS.Vkick[iS]*u.m.to(u.km),
                        SS.costh[iS],np.cos(SS.phi[iS]),SS.Mns[iS]*u.kg.to(u.M_sun),SS.Apost[iS]*u.m.to(u.R_sun),SS.epost[iS],SS.V_sys[iS]*u.m.to(u.km),
                        thatflag,switchflag,disdum,codum,cedum])
    
    return prog 
            
if __name__ == "__main__":
    args=parse_commandline()

    BSE=np.loadtxt(args.BSEoutfile+'.dat')
    temp = BSE.transpose()
    i0 = np.where(temp[9]==0)[0]
    BSE = BSE[i0]
    
    Nb=args.NSNbatch
    size=int(np.floor(len(BSE)/Nb))
    Batches=[BSE[i*size:(i+1)*size] for i in range(int(Nb))]
    #print Batches
    for k,batch in enumerate(Batches):

        if args.ncores>1:
                p=Pool(args.ncores)
                RUNS = p.map(deploy,batch)
                RUNS=filter(partial(is_not,None),RUNS)
                RUNS=np.concatenate(RUNS,axis=1)
        else:
            RUNS=[]
            for b in batch:
            #print b

                temp=SN(b)
                if temp is None:continue
                #print len(temp[0])
                RUNS.append(temp)
                #print np.shape(RUNS)
            #print np.shape(RUNS)
                
            A=[len(np.shape(run)) for run in RUNS]
            RUNS=np.concatenate(RUNS,axis=1)
        RUNS=np.transpose(RUNS)
        np.savez_compressed(args.SNoutfile+'_'+str(k)+'.out',RUNS)
        
