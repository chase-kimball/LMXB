import numpy as np
import matplotlib.pyplot as plt
import Observed
import TrajTools
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
fig=plt.figure(figsize=(20,10))
ax1=fig.add_subplot(1,2,1)
ax2=fig.add_subplot(1,2,2)
x_values=np.arange(-15,15,0.01)
y_values=np.arange(-15,15,0.01)
z_values=np.arange(-15,15,0.01)
X1,Y1=np.meshgrid(x_values,y_values)
pot1=TrajTools.potential(X1,Y1,0)
X2,Z2=np.meshgrid(x_values,z_values)
pot2=TrajTools.potential(X2,0,Z2)
ax1.contourf(X1,Y1,pot1,50,cmap='gray')
ax2.contourf(X2,Z2,pot2,50,cmap='gray')
sysnames=['SCO','LS','Aql','Cyg','XTEJ1118','GROJ1655','LSI61','CEN']
color=['black','orange','blue','red','green','purple','grey','magenta']
i=0
for name in sysnames:
    system=Observed.LMXB_Sys(name)
    ax1.quiver(system.X0,system.Y0,system.U0,system.V0,scale=1500,label=system.fullname,color=color[i])
    ax2.quiver(system.X0,system.Z0,system.U0,system.W0,scale=1500,label=system.fullname,color=color[i])
    i+=1
ax1.axhline(0,color='black',linestyle=':')
ax1.axvline(0,color='black',linestyle=':')
ax1.set_xlim(-15,15)
ax1.set_ylim(-15,15)
ax1.legend()
ax1.set_xlabel('X-Axis of Galactic Plane (kpc)',fontsize=25)
ax1.set_ylabel('Y-Axis of Galactic Plane (kpc)',fontsize=25)
ax2.axhline(0,color='black',linestyle=':')
ax2.axvline(0,color='black',linestyle=':')
ax2.set_xlim(-15,15)
ax2.set_ylim(-15,15)
ax2.legend()
ax2.set_xlabel('X-Axis of Galactic Plane (kpc)',fontsize=25)
ax2.set_ylabel('Z-Axis above Galactic Plane (kpc)',fontsize=25)
position_xy=[-9,1,-6,4.5,-3,6.3,-6,12,-14,3,-8.5,3,-10,5.5,-9,-2.2]
position_xz=[-7.5,-0.5,-4,-2,-3,0.5,-5,0.5,-13.5,1.5,-10.3,0.5,-10.5,-0.5,-7.5,6]
for n in range(0,8):
    system=Observed.LMXB_Sys(sysnames[n])
    ax1.text(position_xy[n*2],position_xy[n*2+1],system.fullname,fontsize=7,color=color[n])
    ax2.text(position_xz[n*2],position_xz[n*2+1],system.fullname,fontsize=7,color=color[n])
plt.savefig('Plots/initial_velocity_potential')