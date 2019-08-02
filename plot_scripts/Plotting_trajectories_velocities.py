# -*- coding: utf-8 -*-
	# Copyright (C) Shrujal Ambati, Richard Ackermann, Chase Kimball (2019)
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
	
	__author__ = ['Shrujal Ambati','Richard Ackermann','Chase Kimball']

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle
import Observed
params = {
        # latex

        # fonts
        'font.family': 'serif',
        #'font.serif': 'Palatino',
        #'font.sans-serif': 'Helvetica',
        #'font.monospace': 'Ubunto Mono',

        # figure and axes
       # 'figure.figsize': (10, 5),
       # 'figure.titlesize': 35,
        'axes.grid': False,
        'axes.titlesize':40,
        #'axes.labelweight': 'bold',
        'axes.labelsize': 30,

        # tick markers
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'xtick.major.size': 10.0,
        'ytick.major.size': 10.0,
        'xtick.minor.size': 3.0,
        'ytick.minor.size': 3.0,

        # legend
        'legend.fontsize': 12,
        'legend.frameon': True,
        'legend.framealpha': 0.5,

        # colors
        'image.cmap': 'viridis',

        # saving figures
        'savefig.dpi': 150
        }

plt.rcParams.update(params)

sys_array = np.array(['XTEJ1118', 'SCO', 'LS', 'Aql', 'Cyg', 'GROJ1655', 'LSI61', 'CEN'])
color=['green','black','orange','yellow','red','purple','0.25','magenta']
for i in range(8):
    plt.figure( figsize=(40,20) )
    test_sys = Observed.LMXB_Sys(sys_array[i])
    R_t_mean, t = test_sys.getTrajectory(1.0, 0.001, test_sys.RK4)
   
    for j in range( 2 ):
        ax = plt.subplot( 1,2 ,j+1)
        with open('Traj_data/' + sys_array[i] + '_R_t', 'rb') as f:
            R_t = pickle.load(f)
        for k in range(len(R_t)):
            ax.plot(R_t[k][3], R_t[k][4+j], color = '0.75')
            ax.set_title(test_sys.fullname)
            ax.set_xlabel('X (kpc)')
            if j == 0:    
                ax.set_ylabel('Y (kpc)')
            elif j == 1:
                ax.set_ylabel('Z (kpc)')
        ax.plot(R_t_mean[3], R_t_mean[4+j], color = 'blue',zorder=500)
        if j==0:
            ax.quiver(test_sys.X0,test_sys.Y0,test_sys.U0,test_sys.V0,color=color[i],scale=1500,zorder=1000)
        if j==1:
            ax.quiver(test_sys.X0,test_sys.Z0,test_sys.U0,test_sys.W0,color=color[i],scale=1500,zorder=1000)
        plt.axhline(0,color='black',linestyle=':')
        plt.axvline(0,color='black',linestyle=':')
    plt.savefig('trajectories_velocity_vectors_'+sys_array[i])
    plt.close()
