#!/usr/bin/env python
# ! ## File: network.py
# ! ## Network class and read edges file.
# ! ## See README.md for more information and use
# !-----------------------------------------------------------------------------
# ! SIS epidemic model algorithm based on the article 
# !           "Optimized Gillespie algorithms for the simulation of Markovian
# !            epidemic processes on large and heterogeneous networks"
# ! Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira
# ! 
# ! Please cite the above cited paper (available at <http://wesleycota.com/ ) as reference
# ! to our code.
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
# ! Date      : 27 Mar 2017
# ! Version   : 1.0
# !-----------------------------------------------------------------------------
# ! See README.md for more details
# ! This code is available at <https://github.com/wcota/dynSIS-py>
# ! For performance, see <https://github.com/wcota/dynSIS> (Fortran implementation)
# ! For NetworkX library, see <https://github.com/wcota/dynSIS-networkx> (NetworkX implementation)

import numpy as np
from tools import *

# Defining a network object
class network(object):
    
    def __init__(self):
        self.size = 0
        self.skk = 0
        self.con = None
        self.ini = None
        self.k = None
        
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

# Read edges from file
def readEdges(fname):
    netw = network()
        
    tmp_con = []
    print_info('Reading file...')
    with open(fname, 'rt') as f:
        for i in f:
            li = list(map(int , i.strip().split(','))) # change ',' if another type
            if li[0] == li[1]:
                print_error('Self-connection found! Verify your data.')
            if li[0] < 1 or li[1] < 1:
                print_error('Vertex id MUST be >= 1. Verify your data.')
            
            netw.size = max(netw.size, li[0], li[1])
            tmp_con.append(li)
            netw.skk += 2
            
    # Now let's build the adjacency list. Numpy is used to use less memory.
    print_info('Building the adjacency list...')
    netw.k = np.zeros(netw.size, np.int)
    netw.ini = np.zeros(netw.size, np.int)
    netw.con = np.zeros(netw.skk, np.int)
    
    pos = 0
    print_info('Calculating degrees...')
    for x, y in tmp_con:
        netw.k[x-1] += 1
        netw.k[y-1] += 1
        
    for i in range(0,netw.size):
        netw.ini[i] = pos
        pos = pos + netw.k[i]
            
    # Now we read again the tmp_con and save data (reindex! Nodes begin at zero)
    print_info('Connecting edges...')
    tmp_pos = np.copy(netw.ini)
    for x,y in tmp_con:
        netw.con[tmp_pos[x-1]] = y - 1
        tmp_pos[x-1] += 1
        netw.con[tmp_pos[y-1]] = x - 1
        tmp_pos[y-1] += 1
        
    return netw
