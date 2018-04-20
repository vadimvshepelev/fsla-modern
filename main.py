import numpy as np

import config as cfg
import problem as problem
import solver.exact as exc
import eos.ideal as eos

print("FSLA hydrocode v.0.1")

U = np.zeros((cfg.N_X+2*cfg.const["N_GHOST_CELLS"], 
              cfg.N_Y+2*cfg.const["N_GHOST_CELLS"], 
			  cfg.N_Z+2*cfg.const["N_GHOST_CELLS"], 
			  cfg.const["CONS_VECT_N_SIZE"]))		  
			  
eos_ideal = eos.EOSIdeal(1.4)
problem_toro_test_1_x = problem.CProblem(eos_ideal, *problem.toro_test_1_x)
exact_solver = exc.CExactRiemannSolver(eos_ideal)

flux = exact_solver.calc_flux(0., 0.)


