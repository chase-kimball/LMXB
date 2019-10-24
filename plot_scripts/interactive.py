import numpy as np
import matplotlib
<<<<<<< HEAD
=======
matplotlib.use('TkAgg')
>>>>>>> 38f287210bce1e1e19a8f36c4c435fc0c662eb09
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
Mhe = 5

ap = np.linspace(20, 70, inputs)
ep = np.linspace(0, .99, inputs)

apx, epy = np.meshgrid(ap, ep)

Apre = apx.flatten()
epre = epy.flatten()

nk = 10
data = []
df = pd.DataFrame()

for i in range(inputs**2):
    print(Mcomp, Mhe, Apre, epre)
    sys = SN.System(Mcomp, Mhe, Apre[i], epre[i], nk)
    sys.SN()
    data.append({'Apre': sys.Apre[0], 'epre': sys.epre[0], 'oldSNflag1': np.sum(sys.oldSNflag1)>=threshold, 
                 'SNflag1': np.sum(sys.SNflag1)>=threshold, 'SNflag2': np.sum(sys.SNflag2)>=threshold,
                 'SNflag3': np.sum(sys.SNflag3)>=threshold, 'SNflag4': np.sum(sys.SNflag4)>=threshold,
                 'oldSNflag4': np.sum(sys.oldSNflag4)>=threshold,
                 'SNflag5': np.sum(sys.SNflag5)>=threshold, 'SNflag6': np.sum(sys.SNflag6)>=threshold, 
                 'SNflag7': np.sum(sys.SNflag7)>=threshold})
    
df = pd.DataFrame(data)

fig = plt.figure()
ax = plt.subplot(1,1,1)

rax = plt.axes([0.05, 0.4, 0.1, 0.15]) #play around til it's readable. These are just coordinates for the box where the buttons are displayed
check = CheckButtons(rax, ('flag1', 'flag2', 'flag3', 'flag4', 'flag5', 'flag6', 'flag7'), [False] * 7) # Name the flags however you want and add as many as you need

onPlot, = ax.plot(df['Apre'], df['epre'],'g.',label='Pass')
test_dictionary = {
                   "flag1": {'on_off': False, 'flag_array': df['SNflag1']},
                   "flag2": {'on_off': False, 'flag_array': df['SNflag2']}, 
                   "flag3": {'on_off': False, 'flag_array': df['SNflag3']},
                   "flag4": {'on_off': False, 'flag_array': df['SNflag4']},
                   "flag5": {'on_off': False, 'flag_array': df['SNflag5']}, 
                   "flag6": {'on_off': False, 'flag_array': df['SNflag6']},
                   "flag7": {"on_off": False, 'flag_array': df['SNflag7']}
                   }

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
    onPlot.set_data(df['Apre'].values[Ion],df['epre'].values[Ion])
    
    plt.draw()
check.on_clicked(func)
plt.show()
