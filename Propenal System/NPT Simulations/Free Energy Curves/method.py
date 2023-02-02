
#-----------------------------------------------------------------------------------------
# packages

import mics as mx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#-----------------------------------------------------------------------------------------
# inputs

conv = 1e-27*24.21726e-3 #conversion factor: atm A^3 --> kcal
Navogadro = 6.02214076e23 # Avogadro's number [atoms/mol]
kb = 1.9872e-3 #Boltzmann constant [kcal/mol/K]
Nat = 512.0 #number of propenal molecules

temp = np.array([120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0]) #sampled temperatures
deltaf_alch = 1.12/4.184/(kb*200.0) #Reduced free energy from PSCP calculations

npoints = 201 #number of states for reweighting
nstates = len(temp) #number of sampled states
deltaG = np.zeros(npoints)

#-----------------------------------------------------------------------------------------
# free energy of the sampled states 

samples_liq = mx.pooledsample()
samples_sol = mx.pooledsample()
mx.verbose = True
for state in range(nstates):

    print(f'state {state}')

    kwargs = {}
    kwargs['T'] = temp[state]
    kwargs['beta'] = 1/(kb*temp[state])  
    
    data = pd.read_csv(f'output_liquid_{int(temp[state])}.out', sep=" ")
    data['U_pot'] = data['PotEng']
    data['PV'] = 1.0*data['Volume']*conv*Navogadro
    data.drop(index=range(5000), inplace=True)
    data.index = np.arange(0, len(data))
    samples_liq.append(mx.sample(data, 'beta*(U_pot + PV)', acfun='U_pot', **kwargs))

    data = pd.read_csv(f'output_solid_{int(temp[state])}.out', sep=" ")
    data['U_pot'] = data['PotEng']
    data['PV'] = 1.0*data['Volume']*conv*Navogadro
    data.drop(index=range(5000), inplace=True)
    data.index = np.arange(0, len(data))
    samples_sol.append(mx.sample(data, 'beta*(U_pot + PV)', acfun='U_pot', **kwargs))

samples_liq.subsampling(integratedACF=True)
samples_sol.subsampling(integratedACF=True)

mixture_liq = mx.mixture(samples_liq, engine=mx.MBAR())
mixture_sol = mx.mixture(samples_sol, engine=mx.MBAR())

results_liq = mixture_liq.free_energies()
results_sol = mixture_sol.free_energies()

#-----------------------------------------------------------------------------------------
# reweighting

variables = dict(beta=[], T=[])
temp_new = np.linspace(temp.min(), temp.max(), npoints)

for point in range(npoints):

    variables['beta'].append( 1/(kb*temp_new[point])  )
    variables['T'].append( temp_new[point]  )

reweighting_liq = mixture_liq.reweighting(
    potential='beta*(U_pot + PV)',
    conditions=pd.DataFrame(variables)
    )

reweighting_sol = mixture_sol.reweighting(
    potential='beta*(U_pot + PV)',
    conditions=pd.DataFrame(variables)
    )

reweighting_liq['f'] = reweighting_liq['f']/Nat
reweighting_sol['f'] = reweighting_sol['f']/Nat
results_liq['f'] = results_liq['f']/Nat
results_sol['f'] = results_sol['f']/Nat

#-----------------------------------------------------------------------------------------
# free energy difference between solid and liquid phases

for point in range(npoints):
    deltaG[point] = (reweighting_sol['f'][point] - reweighting_liq['f'][point] - \
                     results_sol['f'][8] + results_liq['f'][8] + \
                         deltaf_alch)*kb*reweighting_liq['T'][point]

#-----------------------------------------------------------------------------------------
# melting temperature

def integrand_t(x, coefs):
    return coefs[0]*x + coefs[1] 
coefs = np.polyfit(reweighting_liq['T'], deltaG, 1)

root = np.roots(coefs)

#-----------------------------------------------------------------------------------------
# plots

fig = plt.figure(figsize=(3.0, 4.0),  dpi=300)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.set_ylabel(r'$\beta G - \beta_{ref} G_{ref}$')
ax2.set_ylabel(r'$\Delta G_{l,s}$ (kJ/mol)')
ax2.set_xlabel(r'$T$ (K)')

p1, = ax1.plot(reweighting_liq['T'], reweighting_liq['f'], 'c-', label='Liquid')
p2, = ax1.plot(results_liq['T'], results_liq['f'], 'co', label='')
p3, = ax1.plot(reweighting_sol['T'], reweighting_sol['f'], 'm-', label='Solid')
p4, = ax1.plot(results_sol['T'], results_sol['f'], 'mv', label='')

from matplotlib.legend_handler import HandlerTuple
ax1.legend([(p1, p2), (p3, p4)], ['Liquid', 'Solid'], 
               handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2)

ax2.plot(reweighting_liq['T'], deltaG*4.184, 'k-', label='')
ax2.plot(np.linspace(120, 220, 100), np.linspace(0, 0, 100), 'k--', label='')

ax1.set_xticks([120, 140, 160, 180, 200, 220])
ax2.set_xticks([120, 140, 160, 180, 200, 220])
ax1.set_yticks([-20, -15, -10, -5, 0])
ax2.set_yticks([-2, -1, 0, 1, 2])
ax1.set_xticklabels([])

fig.align_ylabels([ax1, ax2])

plt.tight_layout()

plt.savefig('figure_mbar.png', dpi=300)
plt.show()

