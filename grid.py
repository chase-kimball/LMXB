import Observed
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def get_params(system, std_d):
    sys_dict = Observed.getSysDict(system)
    Mbh = [sys_dict['Mbh'][0] + std_d * sys_dict['Mbh'][1], sys_dict['Mbh'][0] - std_d * sys_dict['Mbh'][1]]
    q = [sys_dict['q'][0] + std_d * sys_dict['q'][1], sys_dict['q'][0] - std_d * sys_dict['q'][1]]
    P = sys_dict['P'][0]

    return Mbh, q, P

def make_grid(Mbh_range = [4, 20], q_range = [0, 1], P_range = [0, 40], amount = 10):
    Mbh_RLO = np.linspace(Mbh_range[0], Mbh_range[1], amount)
    q_RLO = np.linspace(q_range[0], q_range[1], amount)
    P_RLO = np.linspace(P_range[0], P_range[1], amount)

    Mbh_RLO, q_RLO, P_RLO = np.meshgrid(Mbh_RLO, q_RLO, P_RLO)

    return Mbh_RLO.flatten(), q_RLO.flatten(), P_RLO.flatten()


def calculate(Mbh_RLO, q_RLO, P_RLO, alpha, beta, params):
    Mdon_RLO = Mbh_RLO * q_RLO

    print('Mbh_RLO:', Mbh_RLO, 'q_RLO:', q_RLO, 'P_RLO:', P_RLO, 'Mdon_RLO:', Mdon_RLO)   

    def equation(v):
        return [(params[2] + P_RLO * (((v[1] * v[0]) / Mdon_RLO) ** (3 * (alpha-1))) * (((Mdon_RLO + Mbh_RLO)/(1+v[1]*v[0]))**2) 
              * math.exp((1-alpha)*(v[1]-(Mdon_RLO/v[0])))), 
               (params[2] + P_RLO * (((v[1] * v[0]) / Mdon_RLO) ** (3 * (alpha-1))) * (((Mdon_RLO + Mbh_RLO)/(1+v[1]*v[0]))**2)
              * ((Mbh_RLO/v[0])**(3/(beta-1))))]

    try: 
        sol = optimize.root(equation, [params[0], params[1]])
        return sol.v
    except:
        return 'error'

main_params = [get_params('GRS1915', 0)[0][0], get_params('GRS1915', 0)[1][0], get_params('GRS1915', 0)[2]]
print(main_params)

grid = make_grid()

for i in range(len(grid[0])):
    print(calculate(grid[0][i], grid[1][i], grid[2][i], 0, 0, main_params))
