import numpy as np
import scipy.linalg as sl
import matplotlib.pyplot as plt
import pandas as pd

import glob, pickle, json

#import hasasia.sensitivity as hsen
#import hasasia.sim as hsim
#import hasasia.skymap as hsky

#import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 300
#mpl.rcParams['figure.figsize'] = [5,3]
#mpl.rcParams['text.usetex'] = True

#from enterprise.pulsar import Pulsar as ePulsar

#import par ant tim files and noise dictionary for the 12y NANOGrav
#pardir = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/12p5yr_stochastic_analysis/data/par/'
#timdir = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/12p5yr_stochastic_analysis/data/tim/'
noise_dir = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/12p5yr_stochastic_analysis/data/'
#pars = sorted(glob.glob(pardir+'*.par'))
#tims = sorted(glob.glob(timdir+'*.tim'))
noise_files = sorted(glob.glob(noise_dir+'*.json'))

#import name of 33 pulsars (note that there are 3 pulsars missing in this list)
data = pd.read_csv('/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/NANOGrav_psr_list')
data_np = data.to_numpy()
data_np_column = data_np[:,0]
psr_list = data_np_column.tolist()

#finalise par and tim files and noise dictionary
#def get_psrname(file,name_sep='_'):
#    return file.split('/')[-1].split(name_sep)[0]

#pars = [f for f in pars if get_psrname(f) in psr_list]
#tims = [f for f in tims if get_psrname(f) in psr_list]
noise = {}
with open(noise_files[2]) as json_file:
    noise = json.load(json_file)

#load pulsar into enterprise.pulsar.Pulsar class instance
#ePsrs = []
#for par,tim in zip(pars,tims):
#    ePsr = ePulsar(par, tim,  ephem='DE436')
#    ePsrs.append(ePsr)
#    print('\rPSR {0} complete'.format(ePsr.name),end='',flush=True)

#define red noise values for 12 years data

key_logA = 'red_noise_log10_A'
key_gamma = 'red_noise_gamma'
log10_A_1 = [value for key, value in noise.items() if key_logA in key]
log10_A_2 = [10**x for x in log10_A_1]
gamma = [value for key, value in noise.items() if key_gamma in key]
name = [key for key, value in noise.items() if key_logA in key]
for i in range(len(name)):
    name[i] = name[i].replace("_red_noise_log10_A", "")
values = [list(t) for t in zip(log10_A_2, gamma)]
rn_psrs = {}
for i in name:
    for j in values:
        rn_psrs[i] = j
        values.remove(j)
        break


print(rn_psrs)

