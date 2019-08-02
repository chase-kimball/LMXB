# -*- coding: utf-8 -*-
# Copyright (C) Shrujal Ambati(2019)
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

__author__ = ['Shrujal Ambati']

import argparse
import Observed
import TrajTools
import numpy as np
import pickle 

parser = argparse.ArgumentParser()
parser.add_argument('--dt', type = float)
parser.add_argument('--T', type = float)
parser.add_argument('--n_runs', type = int)
parser.add_argument('--sys_name', type = str)

args = parser.parse_args()

test_sys = Observed.LMXB_Sys(args.sys_name)

#creates empty lists for Velocity/Position array and time array
Rarr, tarr = [], []

for run in range(args.n_runs):
    test_sys.setRandUVWXYZ()
    R_t, t = test_sys.getTrajectory(args.T, args.dt, test_sys.RK4)

    Rarr.append(R_t)
    tarr.append(t)


with open(args.sys_name + '_R_t', 'wb') as f:
    pickle.dump(Rarr, f)
    
with open(args.sys_name + '_t', 'wb') as f:
    pickle.dump(tarr, f)