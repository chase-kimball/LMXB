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

def make_2d_grid(mdon_range, p_range):
    mdon_all = np.arange(mdon_range[0], mdon_range[1], mdon_range[2])
    p_all = np.arange(p_range[0], p_range[1], p_range[2])
    mdon, p = np.meshgrid(mdon_all, p_all)
    return mdon, p

def calcP(Mbh_RLO, Mdon_RLO, P_RLO, alpha, beta, P_obs):
    P = P_RLO
    Mbh = Mbh_RLO
    q = 0

    if (alpha == 1 or beta == 1):
        Mdon = Mdon_RLO
        Mdon_min = 0
        q = 0
        while True:#(P < P_obs and Mdon > Mdon_min):
            q = Mdon / Mbh
            P = ((P_RLO * (((q * Mbh) / Mdon_RLO) ** (3 * (alpha-1))) * (((Mdon_RLO + Mbh_RLO)/((1+q)*Mbh))**2)) * (math.exp((1-alpha)*(q-(Mdon_RLO/Mbh)))))
            Mdon -= .1
            if P > P_obs:
                return Mdon + .1, q
            elif Mdon < Mdon_min:
                return (0, 0)
    
    else:
        Mbh_max = Mbh_RLO + Mdon_RLO*(1-beta)*(1-alpha)
        while (P < P_obs and Mbh < Mbh_max):
            q = (Mdon_RLO - ((Mbh - Mbh_RLO)/((1-alpha)*(1-beta))))/Mbh
            P = (P_RLO * (((q * Mbh) / Mdon_RLO) ** (3 * (alpha-1))) * (((Mdon_RLO + Mbh_RLO)/((1+q)*Mbh))**2) * (Mbh_RLO/Mbh)**(3/(beta-1)))
            Mbh += .1
        return Mbh - .1, q

def distance(calc, obs, sigma):
    dist = (calc - obs)/sigma
    
    return dist
        

def plot(Mbh_RLO):
    Mdon_RLO, P_RLO = make_2d_grid((0, 16, .2), (0, params[2], .125))
    w = len(P_RLO)
    h = len(P_RLO[0])

    color_mbh = [[0 for r in range(h)] for c in range(w)]
    color_q = [[0 for r in range(h)] for c in range(w)]
    passes = [[0 for r in range(h)] for c in range(w)]

    for i in range(w):
        for x in range(h):
            Mbh, q = calcP(Mbh_RLO, Mdon_RLO[i][x], P_RLO[i][x], alpha, beta, params[2])

            color_mbh[i][x] = (distance(Mbh, params[0][0], params[0][1]))
            color_q[i][x] = (distance(q, params[1][0], params[1][1]))

            if (abs(color_mbh[i][x]) <= 2 and abs(color_q[i][x] <= 2)):
                passes[i][x] = 1
    colors = {'Mbh': color_mbh, 'q': color_q}

    fig, axs = plt.subplots(1, 2)
    fig.set_size_inches(17, 8)
    fig.suptitle('{}\n Mbh_RLO = {}, alpha = {}, beta = {}'.format(system, n, alpha, beta))

    for ax, shade in zip(axs, ('Mbh', 'q')):
        sc = ax.contourf(Mdon_RLO, P_RLO, colors[shade], cmap = "RdBu_r", norm=norm, levels = np.arange(-5, 5.5, .5), extend='both')
        cn = ax.contour(Mdon_RLO, P_RLO, colors[shade], colors='k', levels=[-2, 0, 2], linestyles='solid')
        ax.contour(Mdon_RLO, P_RLO, passes, colors='w', levels = [1], linewidths=2.5)
        ax.clabel(cn, inline=1, fontsize=10, fmt='%1.0f')
        ax.set_title(shade + ' deviation')
        ax.set_xlabel('Mdon_RLO')
        ax.set_ylabel('P_RLO')
        cb = fig.colorbar(sc, ax = ax, ticks = bounds, extend='both')
        cb.set_label('N_sigma')
    
    plt.savefig('div_plots/fig' + str(n), bbox_inches='tight')
    print('plotted ' + str(n))
    plt.close()

system = 'GRS1915'
params = get_params(system)

norm = matplotlib.colors.Normalize(vmin=-5.,vmax=5.)
bounds = range(-5, 6)

alpha = 0
beta = 1

for n in range(2, 17):
    plot(n)
