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

def parse():
    parser = argparse.ArgumentParser()

    parser.add_argument('-x', '--arraynum', type=int)
    parser.add_argument('--inputfile', '-f', type=str, default = '/projects/b1095/samimp/lmxb/work/SN/data/BSEOut_1.h5')
    parser.add_argument('--nkicks', type=int, default = 1000)

    args = parser.parse_args()
    return args

def runSN(arraynum, inputfile, nkicks, rows):

    inputs = pd.read_hdf(inputfile, key='bse')
    cols = list(inputs.columns)

    Mhe = inputs['Mpre']
    Mcomp = inputs['MdonSN']
    Apre = inputs['Apre']
    epre = inputs['epre']

    data = []

    if arraynum is not None:
        rownum = int(len(inputs.index)/20)
        row1 = args.arraynum * rownum 
        lastrow = rownum + row1

        if lastrow > len(inputs.index):
            lastrow = len(inputs.index)
    else:
        row1 = 0
        lastrow = rows

    cols = ['sn_num', 'bin_num', 'M1Zams', 'M2Zams', 'AZams', 'eZams', 'evol_type_preSN', 'tphys', 'kstar_Mpre', 'kstar_MdonSN', 'Mpre', 'MdonSN', 'epre', 'Apre' 'Mns' 'E_ma' 'rpre' 'Vkick', 'costh', 'Apost', 'epost', 'VSx', 'VSy', 'VSz', 'V_sys', 'oldSNflag1', 'SNflag1', 'SNflag2', 'SNflag3', 'SNflag4', 'SNflag5', 'SNflag6', 'SNflag7', 'SNflags']

    for i in range(row1, lastrow):
        sys = SN.System(Mcomp[i], Mhe[i], Apre[i], epre[i], nkicks)
        sys.SN()
        for n in range(nkicks):
            data.append(n, inputs['bin_num'][i], inputs['M1Zams'][i], inputs['M2Zams'][i], inputs['AZams'][i], inputs['eZams'][i], inputs['evol_type_preSN'][i], inputs['tphys'][i], inputs['kstar_Mpre'][i], sys.Mhe[n], sys.Mcomp[n], sys.epre[n], sys.Apre[n], sys.Mns[n], sys.E_ma[n], sys.rpre[n], sys.Vkick[n], sys.costh[n], sys.phi[n], sys.Apost[n], sys.epost[n], sys.VSx[n], sys.VSy, sys.VSz[n], sys.V_sys[n], sys.oldSNflag1[n], sys.SNflag1[n], sys.SNflag2[n], sys.SNflag3[n], sys.SNflag4[n], sys.SNflag5[n], sys.SNflag6[n], sys.SNflag7[n], np.all([item[n] for item in sys.SNflags]))

    SNdata = pd.DataFrame(dict(zip(cols,data)))
    return SNdata
