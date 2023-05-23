import Observed
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize
import pandas as pd

def get_params(system, std_d):
    sys_dict = Observed.getSysDict(system)
    Mbh = [sys_dict['Mbh'][0] + std_d * sys_dict['Mbh'][1], sys_dict['Mbh'][0] - std_d * sys_dict['Mbh'][1]]
    q = [sys_dict['q'][0] + std_d * sys_dict['q'][1], sys_dict['q'][0] - max(0, std_d * sys_dict['q'][1])]
    P = sys_dict['P'][0]

    return Mbh, q, P

def make_grid(Mbh_range = [4, 20], Mdon_range = [0.1, 20], P_range = [0, 40], amount = 10):
    Mbh_RLO = np.linspace(Mbh_range[0], Mbh_range[1], amount)
    Mdon_RLO = np.linspace(Mdon_range[0], Mdon_range[1], amount)
    P_RLO = np.linspace(P_range[0], P_range[1], amount)

    Mbh_RLO, Mdon_RLO, P_RLO = np.meshgrid(Mbh_RLO, Mdon_RLO, P_RLO)

    return Mbh_RLO.flatten(), Mdon_RLO.flatten(), P_RLO.flatten()

def calculate(Mbh_RLO, Mdon_RLO, P_RLO, alpha, beta, mbh_range, q_range, P_obs):
    def equation(Mbh, q):
        return (P_RLO * (((q * Mbh) / Mdon_RLO) ** (3 * (alpha-1))) * (((Mdon_RLO + Mbh_RLO)/((1+q)*Mbh))**2) * (Mbh_RLO/Mbh)**(3/(beta-1))) - P_obs

    ms = np.linspace(mbh_range[0], mbh_range[1], 100)
    qs = np.linspace(q_range[0], q_range[1], 100)
    Ms, Qs = np.meshgrid(ms, qs)
    evaluate = equation(Ms, Qs)
#    print(np.min(evaluate), np.max(evaluate))
    return np.min(evaluate) < 0 and np.max(evaluate) > 0

def check_for_repeats(line, all_data):
    if (line in all_data):
        return False
    else:
        return True

system = 'GRS1915'
params = get_params(system, 2)
grid = make_grid(amount = 15)

data = []
num = 0
real_num = 0

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(len(grid[0])):
    Mbh_RLO = grid[0][i]
    Mdon_RLO = grid[1][i]
    P_RLO = grid[2][i]
    if calculate(Mbh_RLO, Mdon_RLO, P_RLO, 0, 0, params[0], params[1], params[2]):
        num += 1
        new_data = {'System': system, 'Mbh_RLO': Mbh_RLO, 'Mdon_RLO': Mdon_RLO, 'P_RLO': P_RLO}
        if check_for_repeats(new_data, data):
            data.append(new_data)
            real_num += 1

df = pd.DataFrame(data)

print(len(grid[0]))
print(num, real_num)

ax.scatter(df['Mdon_RLO'], df['Mbh_RLO'], df['P_RLO'], color = 'k', alpha=0.5)

ax.set_xlabel('Mdon_RLO')
ax.set_ylabel('Mbh_RLO')
ax.set_zlabel('P_RLO')

plt.show()
