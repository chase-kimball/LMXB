# -*- coding: utf-8 -*-
# Copyright (C) Samuel Imperato (2019)
#
# This file is part of the LMXB package
#
# LMXB is free software: you can redistribute it and/or modify
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
# along with LMXB.  If not, see <http://www.gnu.org/licenses/>.

__author__ = ['Samuel Imperato']



import eccentricPre as SN
import pandas as pd
import numpy as np
import astropy.units as units
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-x', '--arraynum', type=int)
parser.add_argument('--inputfile', '-f', type=str, default = '/projects/b1095/samimp/lmxb/work/SN/data/BSEOut_1.h5')
parser.add_argument('--nkicks', type=int, default = 1000)

args = parser.parse_args()

inputs = pd.read_hdf(args.inputfile, key='bse')
cols = list(inputs.columns)

Mhe = inputs['Mpre']
Mcomp = inputs['MdonSN']
Apre = inputs['Apre']
epre = inputs['epre']

data = []

if args.arraynum is not None:
    rownum = int(len(inputs.index)/20)
    row1 = args.arraynum * rownum 
    lastrow = rownum + row1

    if lastrow > len(inputs.index):
        lastrow = len(inputs.index)
else:
    row1 = 0
    lastrow = len(inputs.index)

for i in range(row1, lastrow):
    sys = SN.System(Mcomp[i], Mhe[i], Apre[i], epre[i], Nkick=args.nkicks)
    sys.SN()
    for n in range(args.nkicks):
        data.append({'sn_num': n, 'bin_num': inputs['bin_num'][i], 'M1Zams': inputs['M1Zams'][i], 'M2Zams': inputs['M2Zams'][i], 'AZams': inputs['AZams'][i], 'eZams': inputs['eZams'][i], 'evol_type_preSN': inputs['evol_type_preSN'][i], 'tphys': inputs['tphys'][i], 'kstar_Mpre': inputs['kstar_Mpre'][i], 'kstar_MdonSN': inputs['kstar_MdonSN'][i], 'Mpre': sys.Mhe[n]*units.kg.to(units.M_sun), 'MdonSN': sys.Mcomp[n]*units.kg.to(units.M_sun), 'epre': sys.epre[n], 'Apre': sys.Apre[n]*units.m.to(units.R_sun), 'Mns': sys.Mns[n]*units.kg.to(units.M_sun), 'E_ma': sys.E_ma[n], 'rpre': sys.rpre[n]*units.m.to(units.R_sun), 'Vkick': sys.Vkick[n]*units.m.to(units.km), 'costh': sys.costh[n], 'phi': sys.phi[n], 'Apost': sys.Apost[n]*units.m.to(units.R_sun), 'epost': sys.epost[n], 'VSx': sys.VSx[n]*units.m.to(units.km), 'VSy': sys.VSy[n]*units.m.to(units.km), 'VSz': sys.VSz[n]*units.m.to(units.km), 'V_sys': sys.V_sys[n]*units.m.to(units.km), 'oldSNflag1': sys.oldSNflag1[n], 'SNflag1': sys.SNflag1[n], 'SNflag2': sys.SNflag2[n], 'SNflag3': sys.SNflag3[n], 'SNflag4': sys.SNflag4[n], 'SNflag5': sys.SNflag5[n], 'SNflag6': sys.SNflag6[n], 'SNflag7': sys.SNflag7[n], 'SNflags': np.all([item[n] for item in sys.SNflags])})
        #if np.all([item[n] for item in sys.SNflags]):
        #    passdata.append(data[-1])
        #else:
        #    faildata.append(data[-1])

df = pd.DataFrame(data)
#passdf = pd.DataFrame(passdata)
#faildf = pd.DataFrame(faildata)

df.to_hdf('./data/SNdata_{}.h5'.format(args.arraynum), key='SN')
#passdf.to_hdf('./data/SNdata_{}.h5'.format(args.arraynum), key='pass')
#faildf.to_hdf('./data/SNdata_{}.h5'.format(args.arraynum), key='fail')
print(df)
