import numpy as np

import config as cfg
import eos.ideal as eos

class CField:
    """Implements full mesh function type: conservative variables vector, setting initial and boundary conditions etc."""
    def __init__(self, problem, NX, NY, NZ):
        self.U = np.zeros((NX+2*cfg.const["N_GHOST_CELLS"], NY+2*cfg.const["N_GHOST_CELLS"], NZ+2*cfg.const["N_GHOST_CELLS"],
                      cfg.const["CONS_VECT_N_SIZE"]))	
        self.NX = NX
        self.NY = NY
        self.NZ = NZ
        dx = (problem.x_max-problem.x_min)/NX
        print(dx)
        self.x_mesh = np.linspace(problem.x_min - cfg.const["N_GHOST_CELLS"]*dx, problem.x_max + cfg.const["N_GHOST_CELLS"]*dx, 
                                  NX+1+2*cfg.const["N_GHOST_CELLS"])
        dy = (problem.y_max-problem.y_min)/NY
        self.y_mesh = np.linspace(problem.y_min - cfg.const["N_GHOST_CELLS"]*dy, problem.y_max + cfg.const["N_GHOST_CELLS"]*dy, 
                                  NY+1+2*cfg.const["N_GHOST_CELLS"])
        dz = (problem.z_max-problem.z_min)/NZ
        self.z_mesh = np.linspace(problem.z_min - cfg.const["N_GHOST_CELLS"]*dy, problem.z_max + cfg.const["N_GHOST_CELLS"]*dz, 
                                  NZ+1+2*cfg.const["N_GHOST_CELLS"])        
    
    def set(self, problem, eos):    	
        """Approximates initial conditions to the mesh function"""
        i_min = cfg.const["N_GHOST_CELLS"]
        j_min = cfg.const["N_GHOST_CELLS"]
        k_min = cfg.const["N_GHOST_CELLS"]
        i_max = i_min + self.NX
        j_max = j_min + self.NY
        k_max = k_min + self.NZ
        ro_l = problem.ro_l
        ro_r = problem.ro_r
        u_l = problem.u_l
        u_r = problem.u_r
        p_l = problem.p_l
        p_r = problem.p_r        
        
        for i in range(i_min, i_max):
            for j in range(j_min, j_max):
                for k in range(k_min, k_max):                    
                    if problem.dir == 'x':                         
                        v_l = 0.
                        w_l = 0.
                        e_l = eos.gete(ro_l, p_l)
                        E_l = e_l + u_l*u_l/2. + v_l*v_l/2. + w_l*w_l/2.
                        v_r = 0.
                        w_r = 0.
                        e_r = eos.gete(ro_r, p_r)
                        E_r = e_r + u_r*u_r/2. + v_r*v_r/2. + w_r*w_r/2.
                        if self.x_mesh[i] < problem.q_0 :
                            self.U[i][j][k] = [ro_l, ro_l*u_l, ro_l*v_l, ro_l*w_l, ro_l*E_l]
                        else:
                            self.U[i][j][k] = [ro_r, ro_r*u_r, ro_r*v_r, ro_r*w_r, ro_r*E_r]
                            
                            
                            
                            
                            