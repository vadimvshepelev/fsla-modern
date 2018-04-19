import numpy as np

import config as cfg
import problem as problem
import solver.exact as exc

print("FSLA hydrocode v.0.1")

U = np.zeros((cfg.N_X+2*cfg.const["N_GHOST_CELLS"], 
              cfg.N_Y+2*cfg.const["N_GHOST_CELLS"], 
			  cfg.N_Z+2*cfg.const["N_GHOST_CELLS"], 
			  cfg.const["CONS_VECT_N_SIZE"]))		  

p = problem.CProblem(*problem.eos_ideal, *problem.toro_test_1_x)
			  
			  
print(exc.calc_flux(cfg.const, 0., 0.))	


