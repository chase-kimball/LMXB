# -*- coding: utf-8 -*-
# Copyright (C) Shrujal Ambati, Richard Ackermann(2019)
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

__author__ = ['Shrujal Ambati','Richard Ackermann']


from matplotlib import *
from matplotlib.patches import ConnectionPatch
import Observed
import TrajTools
import scipy.stats as stats
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy import constants as C
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy import units as u
import numpy as np

gs = gridspec.GridSpec(2,2)
fig = plt.figure(figsize=(20,15))
ax1 = plt.subplot2grid(shape=(2,2), loc=(0,0))
ax2 = plt.subplot2grid(shape=(2,2), loc=(1,0))
ax3 = plt.subplot2grid(shape=(2,2), loc=(1,1))

plt.subplots_adjust(hspace=0)

ax1.axes.get_xaxis().set_visible(False)
ax3.axes.get_yaxis().set_visible(False)
ax3.axes.get_xaxis().set_visible(False)
#gs.tight_layout(h_pad=0.0, w_pad = 0.0)
XTE = Observed.LMXB_Sys('XTEJ1118')
#R_t_XTE,t_XTE=XTE.getTrajectory(1.0, 0.001, XTE.RK4)
R_t_XTE_mean, t_mean = XTE.getTrajectory(1.0, 0.001, XTE.RK4)
with open('XTEJ1118_R_t', 'rb') as f:
    R_t_XTE = pickle.load(f)
for k in range(len(R_t_XTE)):
    ax1.plot(t_mean*u.s.to(u.Gyr), R_t_XTE[k][5], color = 'gray')
    ax2.plot(t_mean*u.s.to(u.Gyr), TrajTools.getPec(R_t_XTE[k][3], R_t_XTE[k][4], R_t_XTE[k][5], R_t_XTE[k][0], R_t_XTE[k][1], R_t_XTE[k][2]), color = 'gray')

ax1.plot(t_mean*u.s.to(u.Gyr), R_t_XTE_mean[5], color = 'blue')
ax2.plot(t_mean*u.s.to(u.Gyr), TrajTools.getPec(R_t_XTE_mean[3], R_t_XTE_mean[4], R_t_XTE_mean[5], R_t_XTE_mean[0], R_t_XTE_mean[1], R_t_XTE_mean[2]), color = 'red')

vPec_XTE = TrajTools.getPec(R_t_XTE_mean[3], R_t_XTE_mean[4], R_t_XTE_mean[5], R_t_XTE_mean[0], R_t_XTE_mean[1], R_t_XTE_mean[2])
i0 = []

for i in range(len(R_t_XTE_mean[5])-1):
    if (R_t_XTE_mean[5][i] > 0 and R_t_XTE_mean[5][i+1] < 0) or (R_t_XTE_mean[5][i] < 0 and R_t_XTE_mean[5][i+1] > 0):
        i0.append(i)
    ax2.scatter(np.array(t_mean*u.s.to(u.Gyr))[i0], np.array(vPec_XTE)[i0], color = 'orange', zorder = 1000)
    
plt.subplots_adjust(wspace=0)
hist_data = np.loadtxt('Trajectories/Vpecs0_XTEJ1118.dat')
ax3.hist(hist_data, orientation = 'horizontal', color='orange', density = True, bins = 75)


ax1.set_xlim(0, min(t_mean*u.s.to(u.Gyr)))
ax2.set_xlim(0, min(t_mean*u.s.to(u.Gyr)))
for j in range(len(i0)):
    coord = (t_mean[i0[j]]*u.s.to(u.Gyr), vPec_XTE[i0[j]])
    con = ConnectionPatch(xyA = coord, xyB = (coord[0],0), coordsA = 'data', coordsB = 'data', axesA = ax2, axesB = ax1, zorder = 10, linewidth = 1.5)
    ax2.add_artist(con)

ax1.set_ylabel('Z (kpc)')
ax2.set_ylabel('Peculiar Velocity (km/s)')
ax2.set_xlabel('time before present (Gyr)')
ax3.set_xlabel('PDF')

ax1.axhline(0,linestyle='--')
plt.savefig('XTE_Lplot')