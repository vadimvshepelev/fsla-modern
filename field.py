# 'field.py' module, CField class implementation

import numpy as np

import config as cfg
import eos.ideal as eos

class CField:
    """Implements full mesh function type: conservative variables vector, setting initial and boundary conditions etc."""
    def __init__(self, problem, eos, NX, NY, NZ):
        print("Class CField: Initializing 3D-mesh function, setting i.c.s and b.c.s...", end="")
        self.U = np.zeros((NX+2*cfg.const["N_GHOST_CELLS"], 
                           NY+2*cfg.const["N_GHOST_CELLS"], 
                           NZ+2*cfg.const["N_GHOST_CELLS"],
                           cfg.const["CONS_VECT_N_SIZE"]))	
        self.U_new = np.copy(self.U)
        self.dx = (problem.x_max-problem.x_min)/NX
        self.x_mesh = np.linspace(problem.x_min - cfg.const["N_GHOST_CELLS"]*self.dx, 
                                  problem.x_max + cfg.const["N_GHOST_CELLS"]*self.dx, 
                                  NX+1+2*cfg.const["N_GHOST_CELLS"])
        self.dy = (problem.y_max-problem.y_min)/NY
        self.y_mesh = np.linspace(problem.y_min - cfg.const["N_GHOST_CELLS"]*self.dy, 
                                  problem.y_max + cfg.const["N_GHOST_CELLS"]*self.dy, 
                                  NY+1+2*cfg.const["N_GHOST_CELLS"])
        self.dz = (problem.z_max-problem.z_min)/NZ
        self.z_mesh = np.linspace(problem.z_min - cfg.const["N_GHOST_CELLS"]*self.dz, 
                                  problem.z_max + cfg.const["N_GHOST_CELLS"]*self.dz, 
                                  NZ+1+2*cfg.const["N_GHOST_CELLS"]) 
        self.i_min = cfg.const["N_GHOST_CELLS"]
        self.j_min = cfg.const["N_GHOST_CELLS"]
        self.k_min = cfg.const["N_GHOST_CELLS"]
        self.i_max = self.i_min + NX
        self.j_max = self.j_min + NY
        self.k_max = self.k_min + NZ		
        self.set_ic(problem, eos)
        self.set_bc(problem)
        print("done!")
    
    def set_ic(self, problem, eos):    	
        """Sets initial conditions to the mesh function everywhere but ghost-cells boundary layers"""
        i_min = self.i_min
        j_min = self.j_min
        k_min = self.k_min
        i_max = self.i_max 
        j_max = self.j_max
        k_max = self.k_max
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
                            
    def set_bc(self, problem):
        """Sets boundary conditions to the 6 boundary layers of config.const['N_GHOST_CELLS'] fictious cells"""
        bcs = problem.bcs
        n_bound = cfg.const['N_GHOST_CELLS']
        # Left X-b.c.
        for i in range(0, self.i_min):
            for j in range(self.j_min, self.j_max):
                for k in range(self.k_min, self.k_max): 
                    if bcs[0] == 't':   
                        self.U[i][j][k] = self.U[self.i_min][j][k]                        
        # Right X-b.c.
        for i in range(self.i_max, self.i_max+n_bound):
            for j in range(self.j_min, self.j_max):
                for k in range(self.k_min, self.k_max): 
                    if bcs[1] == 't':
                        self.U[i][j][k] = self.U[self.i_max-1][j][k]              
        # Left Y-b.c.
        for i in range(0, self.i_max+n_bound):
            for j in range(0, self.j_min):
                for k in range(self.k_min, self.k_max): 
                    if bcs[2] == 't':
                        self.U[i][j][k] = self.U[i][self.j_min][k]                
        # Right Y-b.c.
        for i in range(0, self.i_max+n_bound):
            for j in range(self.j_max, self.j_max+n_bound):
                for k in range(self.k_min, self.k_max): 
                    if bcs[3] == 't':
                        self.U[i][j][k] = self.U[i][self.j_max-1][k]                
        # Left Z-b.c.
        for i in range(0, self.i_max+n_bound):
            for j in range(0, self.j_max+n_bound):
                for k in range(0, self.k_min): 
                    if bcs[4] == 't':
                        self.U[i][j][k] = self.U[i][j][self.k_min]
        # Right Z-b.c.
        for i in range(0, self.i_max+n_bound):
            for j in range(0, self.j_max+n_bound):
                for k in range(0, self.k_max): 
                    if bcs[5] == 't':
                        self.U[i][j][k] = self.U[i][j][self.k_max-1]

    def write_file(self, file_name, t):
        """Dumps mesh function to Tecplot data file for visualization"""
        NX = self.i_max-self.i_min
        NY = self.j_max-self.j_min
        NZ = self.k_max-self.k_min
        print("Function CField.write_file(): writing U field to file " + file_name + " at time t =", t, "...", end="") 
        f = open(cfg.const["OUTPUT_DIR"]+file_name, 'w')
        f.write('VARIABLES="X","Y","Z","ro","ro*u","ro*v","ro*w","ro*E"\n')
        f.write('TITLE="Conservative variables vector field t = ' + str(t) + '"\n')
        f.write('ZONE T="Numerical", I=%d, J=%d, K=%d, F=POINT\n' % (NX, NY, NZ))
        for i in range(self.i_min, self.i_max):
            for j in range(self.j_min, self.j_max):
                for k in range(self.k_min, self.k_max):   
                    f.write("%f %f %f %f %f %f %f %f\n" % (self.x_mesh[i], self.y_mesh[j], self.z_mesh[k], 
                            *self.U[i][j][k]))
        f.close()
        print("done!")
                             
        
                            
                            
                            