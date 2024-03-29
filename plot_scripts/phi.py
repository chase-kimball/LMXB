import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import sys
sys.path.append('..') 
import eccentricPre as SN
import pandas as pd
import numpy as np
import astropy.units as units
inputs = 50

Mcomp = 20 #np.random.choice(a=range(3,21), size=50, replace=false)
Mhe = 5 #np.random.choice(a=range(3,10), size=50, replace=false)
ap = np.linspace(20, 70, inputs)
ep = np.linspace(0.0, 0.5, inputs)

apx, epy = np.meshgrid(ap, ep)

Apre = apx.flatten()
epre = epy.flatten()

nk = 1
data = []
df = pd.DataFrame()

for i in range(inputs**2):
    sys = SN.System(Mcomp, Mhe, Apre[i], epre[i], nk)
    sys.SN()
    for n in range(nk):
        data.append({'phi': sys.phi[n], 'Apre': sys.Apre[n], 'Apost': sys.Apost[n], 'epre': sys.epre[n], 'oldSNflag1': sys.oldSNflag1[n], 'SNflag1': sys.SNflag1[n], 'SNflag2': sys.SNflag2[n], 'SNflag3': sys.SNflag3[n], 'SNflag4': sys.SNflag4[n], 'SNflag5': sys.SNflag5[n], 'SNflag6': sys.SNflag6[n], 'SNflag7': sys.SNflag7[n], 'SNflags': np.all([item[n] for item in sys.SNflags])})
df = pd.DataFrame(data)

fdf = df[df['SNflag4'] == False]
print df['phi']
print fdf['phi']
plt.hist(fdf['phi']) 
plt.title('Phi Where Flag4 Fails')
plt.show()
