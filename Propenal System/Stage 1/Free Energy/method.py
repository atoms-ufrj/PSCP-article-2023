
#-----------------------------------------------------------------------------------------
# packages

import mics as mx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#-----------------------------------------------------------------------------------------
# switching functions and their derivatives

def h1(x):
    return 1.0 if x < 0 else (eta if x > 1 else (1.0-eta)*(-3.0*(x**4)+8.0*(x**3)-6.0*(x**2)) + 1.0)

def h1p(x): #derivative of h1(x)
    return 0.0 if x < 0 else (0.0 if x > 1 else (1.0-eta)*(-12.0*(x**3)+24.0*(x**2)-12.0*(x)))

def h2(x):
    return 1.0 if x < 0 else (eta**2 if x > 1 else (1.0-(eta)**2)*(-3.0*(x**4)+8.0*(x**3)-6.0*(x**2)) + 1.0)

def h2p(x): #derivative of h2(x)
    return 0.0 if x < 0 else (0.0 if x > 1 else (1.0-(eta)**2)*(-12.0*(x**3)+24.0*(x**2)-12.0*(x)))

def h3(x):
    return coeff1*(x**2)*((1.0-x)**3)

def h3p(x): #derivative of h3(x)
    return coeff1*(2*x*((1.0-x)**3) - 3*(x**2)*((1.0-x)**2))

#-----------------------------------------------------------------------------------------
# inputs

T = 200.0 #temperature
kb = 1.9872e-3 #Boltzmann constant [kcal/mol/K]

# switching function coefficients
lbdv = 0.9
lbdc = 0.7
eta = 0.1
coeff1 = 22.0

#sampled lambda states
lambda_v = np.array([0.0001, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.74, 0.75, 0.76, 0.765, 0.775, 0.79, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.9999])

nstates = len(lambda_v) #number of sampled states
npoints = 201 #number of states for reweighting
Nat = 512.0 #number of propenal molecules
kT = T*kb

U = 'h1*U_vdwl + h2*U_coulk + h3*U_gauss' #reduced potential energy
Up = 'h1p*U_vdwl + h2p*U_coulk + h3p*U_gauss' #derivative of the reduced potential energy

#-----------------------------------------------------------------------------------------
# free energy of the sampled states 

dUg = np.zeros(nstates)
dUv = np.zeros(nstates)
dUc = np.zeros(nstates)
samples = mx.pooledsample()
mx.verbose = True

for state in range(nstates):
    print(f'state {state}')

    lbda = lambda_v[state]

    data = pd.read_csv(f'output_{state}.out', sep=" ")


    data['U_vdwl'] = data['c_lj']/h1(lbda/lbdv) + data['E_tail']/h1(lbda/lbdv)
    data['U_coulk'] = data['c_coul']/h2(lbda/lbdc) + data['E_long']/h2(lbda/lbdc) 
    data['U_gauss'] = data['c_gauss']/h3(lbda)
    
    data['dUg'] = data['U_gauss']*h3p(lbda)
    data['dUvdwl'] = data['U_vdwl']*h1p(lbda/lbdv)
    data['dUcoulk'] = data['U_coulk']*h2p(lbda/lbdc)

    data.drop(index=range(5000), inplace=True)
    data.index = np.arange(0, len(data))

    kwargs = {'lbda': lbda, 'beta': 1.0/kT}
    kwargs['h1'] = h1(lbda/lbdv)   
    kwargs['h2'] = h2(lbda/lbdc) 
    kwargs['h3'] = h3(lbda) 
 
    samples.append(mx.sample(data, f'beta*({U})', acfun='PotEng+KinEng', **kwargs))
    
    dUv[state] = np.mean(data['dUvdwl'].values)
    dUc[state] = np.mean(data['dUcoulk'].values)
    dUg[state] = np.mean(data['dUg'].values)
    
samples.subsampling(integratedACF=True)

mixture = mx.mixture(samples, engine=mx.MBAR())

results = mixture.free_energies()
results['F'] = results['f']*kT
results['dF'] = results['df']*kT

#-----------------------------------------------------------------------------------------
# reweighting

variables = dict(lbda=[], h1=[], h2=[], h3=[], h1p=[], h2p=[], h3p=[])
for point in range(npoints):

    lbda = point/(npoints - 1)

    variables['lbda'].append(lbda)

    variables['h1'].append(h1(lbda/lbdv))
    variables['h1p'].append(h1p(lbda/lbdv))
    variables['h2'].append(h2(lbda/lbdc))
    variables['h2p'].append(h2p(lbda/lbdc))
    variables['h3'].append(h3(lbda))
    variables['h3p'].append(h3p(lbda))

properties = dict(
    Up=f'{Up}'
    )

combinations = dict(
    F='kT*f'   
    )

reweighting = mixture.reweighting(
    potential=f'beta*({U})',
    properties=properties,
    combinations=combinations,
    conditions=pd.DataFrame(variables), 
    beta=1/kT, 
    kT=kT
    )

#-----------------------------------------------------------------------------------------
# final free energy

deltaA = results['F'].iloc[-1]/Nat*4.184
d_deltaA = results['dF'].iloc[-1]/Nat*4.184

print(f'Delta A = {deltaA} {d_deltaA} kJ/mol')

#-----------------------------------------------------------------------------------------
# plots

fig, ax = plt.subplots(2, 1, figsize=(2.5, 5.0),  dpi=300)
ax[0].plot(reweighting['lbda'], reweighting['F']/Nat*4.184, 'k-', label='Interpolated')
ax[0].plot(results['lbda'], results['F']/Nat*4.184, 'bo', label='Simulated')
ax[0].legend()
ax[1].plot(reweighting['lbda'], reweighting['Up']/Nat*4.184, 'k-', label='')
ax[1].plot(lambda_v, dUv/Nat*4.184 + dUc/Nat*4.184 + dUg/Nat*4.184, 'bo', label='')
ax[0].set_ylabel('$A$ (kJ/mol)')
ax[1].set_ylabel(r'$\partial A/\partial \lambda$ (kJ/mol)')
ax[1].set_xlabel('$\lambda$')
ax[0].set_xticklabels([])
plt.savefig('figure_MBAR_total.png', dpi=300)
plt.show()

