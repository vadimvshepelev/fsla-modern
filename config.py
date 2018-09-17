# CONS_VECT_N_SIZE -- number of elements in conservative variables vector. Typically equals 5 for one-component 3D problems
# N_GHOST_CELLS -- number of boundary ghost cells, depends on the approximation order of the method. 1 for 1st-order, 2 for 2nd-order etc. 

const = {"CONS_VECT_N_SIZE": 5,
         "N_GHOST_CELLS": 2,
         "OUTPUT_DIR": "./output/",
         "CFL": .3,
         "NX": 200,
         "NY": 20,
         "NZ": 20}

