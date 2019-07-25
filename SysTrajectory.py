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