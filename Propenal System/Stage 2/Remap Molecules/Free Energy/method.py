
#-----------------------------------------------------------------------------------------
# packages

import mics as mx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#-----------------------------------------------------------------------------------------
# inputs

T = 200.0 #temperature [K]
kb = 1.9872e-3 #Boltzmann constant [kcal/mol/K]
conv = 1e-27*24.21726e-3 #atm A3 --> kcal
Navogadro = 6.02214076e23 #[atoms/mol]
Nat = 512.0 #number of propenal molecules
nstates = 21 #number of sampled states
kT = T*kb

#-----------------------------------------------------------------------------------------
# free energy of the sampled states 

Volume = np.zeros(nstates)
samples = mx.pooledsample()
mx.verbose = True

for state in range(nstates):
    print(f'state {state}')
     
    data = pd.read_csv(f'State-{state}/{state}_output_{state}.txt', sep=" ")  
    
    for i in range(nstates):
        data_new = pd.read_csv(f'State-{state}/{state}_output_{i}.txt', sep=" ") 
        data[f'Upot_{i}'] = data_new['PotEng'] 
        
    data.index = np.arange(0, len(data))
    Volume[state] = data['Volume'][0]
    kwargs = {'beta': 1.0/kT}
    kwargs['N'] = 512.0
    kwargs['V'] = Volume[state]
    samples.append(mx.sample(data, f'-N*ln(V)+beta*(Upot_{state})', acfun='PotEng', **kwargs)) 
    
samples.subsampling(integratedACF=True)

mixture = mx.mixture(samples, engine=mx.MBAR())

results = mixture.free_energies()
results['F'] = results['f']*kT
results['dF'] = results['df']*kT

#-----------------------------------------------------------------------------------------
# final free energy

deltaA = results['F'].iloc[-1]/Nat*4.184
d_deltaA = results['dF'].iloc[-1]/Nat*4.184

print(f'Delta A = {deltaA} {d_deltaA} kJ/mol')

#-----------------------------------------------------------------------------------------
# plot

fig = plt.figure(figsize=(2.5, 2.5),  dpi=300)
plt.plot(results['V'], results['F']/Nat*4.184, 'k-', label='')
plt.plot(results['V'], results['F']/Nat*4.184, 'bo', label='')
plt.ylabel('$A$ (kJ/mol)')
plt.xlabel('$V$ (\AA$^3$)')
plt.savefig('figure_MBAR_total.png')
plt.show()


