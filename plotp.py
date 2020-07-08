import numpy as np
import Observed
import matplotlib.pyplot as plt
import math
def calculate_P(alpha, beta, P_RLO, Mbh_RLO, Mdon_RLO, amount):
    P = []
    Mbh_max = Mbh_RLO + Mdon_RLO*(1-beta)*(1-alpha)
    Mbh = np.linspace(Mbh_RLO, Mbh_max, amount)

    print(Mbh_max)

    for i in range(amount):
        q = (Mdon_RLO - ((Mbh[i] - Mbh_RLO)/((1-alpha)*(1-beta))))/Mbh[i]
        P.append(P_RLO * (((q * Mbh[i]) / Mdon_RLO) ** (3 * (alpha-1))) * (((Mdon_RLO + Mbh_RLO)/((1+q)*Mbh[i]))**2) * (Mbh_RLO/Mbh[i])**(3/(beta-1)))
    return P, Mbh

amount = 100

P_obs = Observed.getSysDict('GRS1915')['P'][0]
Mbh_obs = Observed.getSysDict('GRS1915')['Mbh'][0]

#P, Mbh = calculate_P(0, 0, 3, 9, 0.5, amount)
#plt.plot(Mbh, P)

for b in [0.8, 0.5, 0]:
    P, Mbh = calculate_P(0, b, 3, 9, 0.5, amount)
    plt.plot(Mbh, P, label='beta = ' + str(b))

plt.plot(Mbh, [P_obs]*amount, label='P observed')

plt.xlabel('Mbh')
plt.ylabel('P')
plt.ylim(0, 100)
plt.legend()
plt.show()
