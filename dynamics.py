#!/usr/bin/env python
# ! ## File: dynamics.py
# ! ## See README.md for more information and use
# !-----------------------------------------------------------------------------
# ! SIS epidemic model algorithm based on the article "Simulation of Markovian epidemic models on large networks"
# ! Copyright (C) 2016 Wesley Cota, Silvio C. Ferreira
# ! 
# ! Please cite the above cited paper as reference to our code.
# ! 
# !    This program is free software: you can redistribute it and/or modify
# !    it under the terms of the GNU General Public License as published by
# !    the Free Software Foundation, either version 3 of the License, or
# !    (at your option) any later version.
# !
# !    This program is distributed in the hope that it will be useful,
# !    but WITHOUT ANY WARRANTY; without even the implied warranty of
# !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# !    GNU General Public License for more details.
# !
# !    You should have received a copy of the GNU General Public License
# !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# !-----------------------------------------------------------------------------
# ! Author    : Wesley Cota
# ! Email     : wesley.cota@ufv.br
# ! Date      : October 2016
# !-----------------------------------------------------------------------------
# ! See README.md for more details
# ! This code is available at <https://github.com/wcota/dynSIS-py>

from network import *
from tools import *
import sys

print_header()

# READING PARAMETERS
if len(sys.argv) < 2:
    print_error('You must enter input and output names as arguments!')

fnInput = sys.argv[1]
fnOutput = sys.argv[2]

print_info('Reading dynamical parameters...')
    
dynp_sam = int(input('How much dynamics samples? '))
dynp_lb = float(input('Value of infection rate lambda (mu is defined as equal to 1) '))
dynp_tmax = int(input('Maximum time steps (it stops if the absorbing state is reached) '))
dynp_pINI = float(input('Fraction of infected vertices on the network as initial condition (is random to each sample) '))
# / READING PARAMETERS

# LOADING NETWORK
print_info('Loading network to memory...')

netw = readEdges(fnInput)

print_info('Everything ok!')
# / LOADING NETWORK

# PREPARING THE NECESSARY THINGS
net_kmax = max(netw.k)
avg_rho = numpy.zeros(dynp_tmax, numpy.float64)
avg_t = numpy.zeros(dynp_tmax, numpy.float64)
avg_sam = numpy.zeros(dynp_tmax, numpy.int)

dyn_ocp = numpy.zeros(netw.size, numpy.int)
dyn_sig = numpy.zeros(netw.size, numpy.int)
# / PREPARING THE NECESSARY THINGS

# RUNNING DYNAMICS
print_info('Running dynamics...',True)

dyn_dt_pos_max = 0 # Auxiliar
for sam in range(1,dynp_sam+1):
    print_info('Sample #'+str(sam),True)

    # Initial conditions
    print_info('Initial condition...')
    dyn_sig[:] = 0.0
    dyn_ocp[:] = 0.0
    dyn_voc = 0
    dyn_sk = 0
    for i in range(0, int(netw.size*dynp_pINI)):
        isOk = False
        while not isOk:
            ver = numpy.random.randint(0,netw.size)
            if dyn_sig[ver] == 0:
                dyn_ocp[dyn_voc] = ver
                dyn_voc += 1
                dyn_sig[ver] = 1
                dyn_sk += netw.k[ver]
                isOk = True
    
    # RUN, Forest, RUN!
    dyn_t = 0
    dyn_dt = 0.0
    dyn_dt_pos = 1
    
    print_info('Running...')
    # SIS II ALGORITHM
    while dyn_t <= dynp_tmax and dyn_voc > 0:
        dyn_p = 1.0*dyn_voc / (dyn_voc + 1.0*dynp_lb * dyn_sk)
        rnd = numpy.random.uniform()
        
        if rnd < dyn_p: # Cure
            # Select a random occupied vertex and heal.
            pos_ocp = numpy.random.randint(0,dyn_voc)
            ver = dyn_ocp[pos_ocp]
                
            # Healed
            dyn_sig[ver] = 0
            dyn_sk -= netw.k[ver]
            dyn_voc -= 1
            dyn_ocp[pos_ocp] = dyn_ocp[dyn_voc]
        else: # Try to infect
            # Select an infected vertex and accept with the probability.
            isOk = False
            while not isOk:
                pos_ocp = numpy.random.randint(0,dyn_voc)
                ver = dyn_ocp[pos_ocp]
                if numpy.random.uniform() < 1.0*netw.k[ver] / (1.0*net_kmax):
                    isOk = True
            
            # Select one of its neighbors
            pos_nei = numpy.random.randint(netw.ini[ver], netw.ini[ver] + netw.k[ver])
            ver = netw.con[pos_nei]
            
            if dyn_sig[ver] == 0:
                dyn_sig[ver] = 1
                dyn_sk += netw.k[ver]
                dyn_ocp[dyn_voc] = ver
                dyn_voc += 1
    
        if dyn_voc > 0:
            dyn_dt = 1.0/(dyn_voc + 1.0*dynp_lb * (1.0*dyn_sk))
            dyn_t += dyn_dt
            
            while (dyn_t >= dyn_dt_pos):
                dyn_dt_pos_max = max(dyn_dt_pos,dyn_dt_pos_max)
                avg_rho[dyn_dt_pos - 1] += 1.0*dyn_voc/netw.size
                avg_t[dyn_dt_pos - 1] += dyn_t
                avg_sam[dyn_dt_pos - 1] += 1
                dyn_dt_pos += 1
    
    # Write output file
    flOutput = open(fnOutput, 'wt')
    print(  '#@ Network file: '+fnInput,
            '#@ Number of nodes: '+str(netw.size),
            '#@ Number of edges: '+str(netw.skk),
            '#@ Samples: '+str(dynp_sam),
            '#! Infection rate (lambda): '+str(dynp_lb),
            '#! Maximum time steps: '+str(dynp_tmax),
            '#! Fraction of infected vertices (initial condition): '+str(dynp_pINI),
            sep='\n',
            file=flOutput)
            
    for dt_pos in range(0,dyn_dt_pos_max):
        print(1.0*avg_t[dt_pos]/avg_sam[dt_pos], 1.0*avg_rho[dt_pos]/(1.0*sam),
                file=flOutput)
                
    flOutput.close()
# / RUNNING DYNAMICS

print_info('Everything ok!',True)
print_info('Input file (edges list): '+ fnInput)
print_info('Output file: '+ fnOutput)