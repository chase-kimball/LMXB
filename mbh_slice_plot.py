import Observed
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn

def get_params(system):
    sys_dict = Observed.getSysDict(system)
    Mbh = (sys_dict['Mbh'][0], sys_dict['Mbh'][1])
    q = (sys_dict['q'][0], sys_dict['q'][1])
    P = sys_dict['P'][0]

    return Mbh, q, P

def make_2d_grid(p1_range, p2_range):
    param1 = np.arange(p1_range[0], p1_range[1], p1_range[2])
    param2 = np.arange(p2_range[0], p2_range[1], p2_range[2])
    p1, p2 = np.meshgrid(param1, param2)
    return p1, p2

def calcP(Mbh_RLO, Mdon_RLO, P_RLO, alpha, beta, P_obs):
    Mbh_max = Mbh_RLO + Mdon_RLO*(1-beta)*(1-alpha)
    P = P_RLO
    Mbh = Mbh_RLO
    q = 0
    
    while (P < P_obs and Mbh < Mbh_max):
        q = (Mdon_RLO - ((Mbh - Mbh_RLO)/((1-alpha)*(1-beta))))/Mbh
        P = (P_RLO * (((q * Mbh) / Mdon_RLO) ** (3 * (alpha-1))) * (((Mdon_RLO + Mbh_RLO)/((1+q)*Mbh))**2) * (Mbh_RLO/Mbh)**(3/(beta-1)))
        Mbh += .1
    return Mbh - .1, q

def distance(calc, obs, sigma):
    dist = (calc - obs)/sigma
    
    return dist
        
system = 'V404'
params = get_params(system)

norm = matplotlib.colors.Normalize(vmin=-5.,vmax=5.)
bounds = range(-5, 6)
#bounds = [-4, -3, -2, -1, 0, 1, 2, 3, 4]


for n in range(2, 17):
    Mbh_RLO = n

    grid = make_2d_grid((0, 16, .2), (0, params[2], .125))
    Mdon_RLO = grid[0]
    P_RLO = grid[1]
    
    w = len(Mdon_RLO)
    h = len(Mdon_RLO[0])

    color_mbh = [[0 for r in range(h)] for c in range(w)] 
    color_q = [[0 for r in range(h)] for c in range(w)] 
    passes = [[0 for r in range(h)] for c in range(w)]

    for i in range(w):
        for x in range(h):
            Mbh, q = calcP(Mbh_RLO, Mdon_RLO[i][x], P_RLO[i][x], 0, 0, params[2])
        
            color_mbh[i][x] = (distance(Mbh, params[0][0], params[0][1]))
            color_q[i][x] = (distance(q, params[1][0], params[1][1]))

            if (abs(color_mbh[i][x]) <= 2 and abs(color_q[i][x]) <= 2):
                passes[i][x] = 1

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(17, 8)
    fig.suptitle('Mbh_RLO: ' + str(n))

    sc1 = ax1.pcolormesh(Mdon_RLO, P_RLO, color_mbh, cmap="RdBu_r", norm=norm)
    cn1 = ax1.contour(Mdon_RLO, P_RLO, color_mbh, colors='k', levels=(-2, 0, 2), linestyles='solid')
    ax1.contour(Mdon_RLO, P_RLO, passes, colors='w', levels=[1], linestyles='dashed')
    ax1.clabel(cn1, inline=1, fontsize=10, fmt='%1.0f')
    ax1.set_title('Mbh deviation')
    ax1.set_xlabel('Mdon_RLO')
    ax1.set_ylabel('P_RLO')
    cb1 = fig.colorbar(sc1, ax = ax1, ticks = bounds, extend='both')
    cb1.set_label('standard deviations away')

    sc2 = ax2.pcolormesh(Mdon_RLO, P_RLO, color_q, cmap="RdBu_r", norm=norm)
    cn2 = ax2.contour(Mdon_RLO, P_RLO, color_q, colors='k', levels=(-2, 0, 2), linestyles='solid')
    ax2.contour(Mdon_RLO, P_RLO, passes, colors='w', levels=[1], linestyles='dashed')
    ax2.clabel(cn2, inline=1, fontsize=10, fmt='%1.0f')
    ax2.set_title('q deviation')
    ax2.set_xlabel('Mdon_RLO')
    ax2.set_ylabel('P_RLO')
    cb2 = fig.colorbar(sc2, ax = ax2, ticks = bounds, extend='both')
    cb2.set_label('standard deviations away')
    
    plt.savefig('div_plots/fig' + str(n), bbox_inches='tight')
    print('plotted ' + str(n))
    plt.close()
