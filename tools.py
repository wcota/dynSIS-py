#!/usr/bin/env python
# ! ## File: tools.py
# ! ## Print and read functions.
# ! ## See README.md for more information and use
# !-----------------------------------------------------------------------------
# ! SIS epidemic model algorithm based on the article
# !           Computer Physics Communications 219C (2017) pp. 303-312
# !           "Optimized Gillespie algorithms for the simulation of 
# !            Markovian epidemic processes on large and heterogeneous networks"
# ! Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira
# ! 
# ! Please cite the above cited paper (available at <http://dx.doi.org/10.1016/j.cpc.2017.06.007> ) 
# ! as reference to our code.
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

# Print functions
def print_error(st):
    raise ValueError(st)

def print_info(st,nl=False):
    if nl:
        print('')
    print('$!', st)
    