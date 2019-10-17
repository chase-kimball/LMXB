import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
import eccentricPre as SN
import astropy.units as units
inputs = 50

Mcomp = 20
Mhe = 5
ap = np.linspace(20, 70, inputs)
ep = np.linspace(20, 70, inputs)

apx, epy = np.meshgrid(ap, ep)

Apre = apx.flatten()
epre = epy.flatten()

nk = 10
data = []
df = pd.DataFrame()

for i in range(inputs**2):
    sys = SN.System(Mcomp, Mhe, Apre[i], epre[i], nk)
    sys.SN()
    data.append({'Apre': sys.Apre, 'epre': sys.epre, 'oldSNflag1': sys.oldSNflag1, 'SNflag1': sys.SNflag1, 'SNflag2': sys.SNflag2, 'SNflag3': sys.SNflag3, 'SNflag4': sys.SNflag4, 'SNflag5': sys.SNflag5, 'SNflag6': sys.SNflag6, 'SNflag7': sys.SNflag7, 'SNflags': np.all([item for item in sys.SNflags])})
df = pd.DataFrame(data)

fig = plt.figure()
ax = plt.subplot(1,1,1)

rax = plt.axes([0.05, 0.4, 0.1, 0.15]) #play around til it's readable. These are just coordinates for the box where the buttons are displayed
check = CheckButtons(rax, ('flag1', 'flag2', 'flag3', 'flag4', 'flag5', 'flag6', 'flag7'), (False) * 7) # Name the flags however you want and add as many as you need

onPlot, = ax.scatter(df['Apre'], df['epre'],'g.',label='Pass')
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
    
    test_dictionary[label]['on_off'] = not test_dictionary[label]['on_off']
    
    Ion = np.ones_like(df['Apre']) #(or like any of your arrays)
    for key in test_dictionary.keys():
        if test_dictionary[key]['on_off']:
            Ion *= test_dictionary[key]['flag_array']

    Ion = np.where(Ion)

    onPlot.set_data(df['Apre'][Ion],df['Apre'][Ion])
    
    plt.draw()
check.on_clicked(func)    
