import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
import astropy.units as units
import sys
sys.path.append('..')
import eccentricPre as SN
import pandas as pd
inputs = 50
threshold = 1
Mcomp = 20
Mhe = 15

ap = np.linspace(20, 70, inputs)
ep = np.linspace(0, .99, inputs)

apx, epy = np.meshgrid(ap, ep)

Apre = apx.flatten()
epre = epy.flatten()

nk = 1
data = []
df = pd.DataFrame()

flags = ['oldSNflag1', 'SNflag1', 'SNflag2', 'SNflag3', 'SNflag4', 'SNflag5', 'SNflag6', 'SNflag7', 'ApostGT0']

for i in range(inputs**2):
    print(Mcomp, Mhe, Apre, epre)
    sys = SN.System(Mcomp, Mhe, Apre[i], epre[i], nk)
    sys.SN()
    data.append({'Apre': sys.Apre[0], 'Apost': sys.Apost[0], 'epre': sys.epre[0],
                 flags[0]: np.sum(sys.oldSNflag1)>=threshold, 
                 flags[1]: np.sum(sys.SNflag1)>=threshold, flags[2]: np.sum(sys.SNflag2)>=threshold,
                 flags[3]: np.sum(sys.SNflag3)>=threshold, flags[4]: np.sum(sys.SNflag4)>=threshold,
                 flags[5]: np.sum(sys.SNflag5)>=threshold, flags[6]: np.sum(sys.SNflag6)>=threshold, 
                 flags[7]: np.sum(sys.SNflag7)>=threshold, flags[8]: sys.Apost[0] > 0})
    
df = pd.DataFrame(data)

fig = plt.figure()
ax = plt.subplot(1,1,1)

rax = plt.axes([0.01, 0.4, 0.1, 0.15]) #play around til it's readable. These are just coordinates for the box where the buttons are displayed
check = CheckButtons(rax, flags, [False] * 9) # Name the flags however you want and add as many as you need

Plot, = ax.plot(df['Apre']/df['Apost'], df['epre'],'k.',label='Fail')
onPlot, = ax.plot(df['Apre']/df['Apost'], df['epre'],'g.',label='Pass')

test_dictionary = {key: {'on_off': False, 'flag_array': df[key]} for key in flags}
 
def func(label):
    global Ion
    global Ioff
    print(label)
    print('before')
    for key in test_dictionary.keys():
        print (key, test_dictionary[key]['on_off'])
    test_dictionary[label]['on_off'] = not test_dictionary[label]['on_off']
    print('after')
    for key in test_dictionary.keys():
        print (key, test_dictionary[key]['on_off'])
    Ion = np.ones_like(df['Apre']) #(or like any of your arrays)

    for key in test_dictionary.keys():
        if test_dictionary[key]['on_off']:
            print (key)
            Ion *= test_dictionary[key]['flag_array']
    Ion = np.where(Ion)[0]
    print(df['Apre'].values[Ion])
    onPlot.set_data(df['Apre'].values[Ion]/df['Apost'].values[Ion],df['epre'].values[Ion])
 
    plt.draw()
check.on_clicked(func)
plt.xlabel('ApreSN/ApostSN')
plt.ylabel('EpreSN')
plt.show()
