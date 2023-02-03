import os, sys, subprocess
import numpy as np

restart = 2000200

ai = np.array([29.099537, 29.152948, 29.206359, 29.25977, 29.31318, 29.366591, 29.420002, 29.473413, 29.526824, 29.580235, 29.633646, 29.687057, 29.740468, 29.793879, 29.84729, 29.900701, 29.954112, 30.007523, 30.060934, 30.114345, 30.167756])
bi = np.array([37.18353, 37.234978, 37.286426, 37.337874, 37.389322, 37.44077, 37.492218, 37.543666, 37.595114, 37.646562, 37.69801, 37.749458, 37.800906, 37.852354, 37.903803, 37.955251, 38.006699, 38.058147, 38.109595, 38.161043, 38.212491])
ci = np.array([42.955202, 43.019741, 43.084281, 43.148821, 43.213361, 43.277901, 43.342442, 43.406983, 43.471524, 43.536065, 43.600607, 43.665149, 43.729692, 43.794235, 43.858778, 43.923321, 43.987865, 44.052408, 44.116953, 44.181497, 44.246042])
alphai = np.array([89.9537, 89.956082, 89.958457, 89.960825, 89.963185, 89.965538, 89.967884, 89.970223, 89.972555, 89.97488, 89.977198, 89.979509, 89.981813, 89.984111, 89.986401, 89.988684, 89.990961, 89.993231, 89.995494, 89.99775, 90.0])
betai = np.array([90.0823, 90.078068, 90.073848, 90.069641, 90.065447, 90.061265, 90.057095, 90.052938, 90.048793, 90.044661, 90.040541, 90.036433, 90.032337, 90.028253, 90.024181, 90.020121, 90.016074, 90.012038, 90.008013, 90.004001, 90.0])
gammai = np.array([90.0155, 90.014705, 90.013912, 90.013121, 90.012332, 90.011545, 90.010761, 90.009978, 90.009198, 90.00842, 90.007644, 90.00687, 90.006099, 90.005329, 90.004562, 90.003796, 90.003033, 90.002272, 90.001512, 90.000755, 90.0])


for j in range(5000):
    for i in range(21): 

        template ='''
variable        NAME index out
log             ${NAME}.log

variable	 a equal $variavel1
variable	 b equal $variavel2
variable	 c equal $variavel3
variable	 alpha equal $variavel4
variable	 beta equal $variavel5
variable	 gamma equal $variavel6

variable        xv equal v_a
variable        xyv equal v_b*cos(v_gamma*PI/180)
variable        yv equal sqrt(v_b*v_b-v_xyv*v_xyv)
variable        xzv equal v_c*cos(v_beta*PI/180)
variable        yzv equal (v_b*v_c*cos(v_alpha*PI/180)-v_xyv*v_xzv)/v_yv
variable        zv equal sqrt(v_c*v_c-v_xzv*v_xzv-v_yzv*v_yzv)

units           real
special_bonds   lj 0.0 0.0 0.5 coul 0.0 0.0 0.8333
#neighbor        2.0 bin
#neigh_modify    delay 0 every 1 check yes page 1000000 one 20000 
atom_style      full
pair_style	 hybrid/overlay lj/cut 12.0 coul/long 12.0 

bond_style	 harmonic
dihedral_style	 charmm
angle_style	 harmonic
improper_style	 cvff

read_restart    restart.out.2.$variavel7
''' 
        template_mod = template*1
        template_mod = template_mod.replace('$variavel1',f'{ai[i]}')
        template_mod = template_mod.replace('$variavel2',f'{bi[i]}')
        template_mod = template_mod.replace('$variavel3',f'{ci[i]}')
        template_mod = template_mod.replace('$variavel4',f'{alphai[i]}')
        template_mod = template_mod.replace('$variavel5',f'{betai[i]}')
        template_mod = template_mod.replace('$variavel6',f'{gammai[i]}')
        template_mod = template_mod.replace('$variavel7',f'{int(restart) + j*200}')

        f=open('in.input','w')
        f.write(template_mod)
        f.close()

        msg  = subprocess.run(['mpirun', '-np', '24', '/afs/crc.nd.edu/user/g/gcorrea2/lammps-29Sep2021/src/lmp_mpi', '-in', 'in.step2.pos'],capture_output=True)
        #msg  = subprocess.run(['/home/gabriela/lammps-29Sep2021/src/lmp_serial', '-i', 'in.step3.pos'],capture_output=True)
        #print(msg.stdout.decode())
        #print(msg.stderr.decode())

        f = open("output.out", "w")
        msg  = subprocess.run(['/afs/crc.nd.edu/user/g/gcorrea2/postlammps/postlammps', '-in', 'ar.log', 'print', 'Step', 'Volume', 'Press', 'KinEng', 'PotEng', 'Enthalpy', 'E_mol', 'E_vdwl', 'E_tail', 'E_coul', 'E_long'], stdout=f)
        #msg  = subprocess.run(['postlammps', '-in', 'ar.log', 'print', 'Step', 'Volume', 'Press', 'KinEng', 'PotEng', 'Enthalpy', 'E_mol', 'E_vdwl', 'E_tail', 'E_coul', 'E_long'], stdout=f)
        f.close()

        fi=open('output.out','r') 
        ft=open(f'20_output_{i}.txt','a') 
    
        fi.readline() #descartando primeira linha
        ft.write(fi.readline())
        fi.close()
        ft.close()
    
    
