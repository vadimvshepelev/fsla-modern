# 'field.py' module, CField class implementation

import numpy as np
import math

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
        self.F = np.copy(self.U)
        self.G = np.copy(self.U)
        self.H = np.copy(self.U)
        self.S = np.copy(self.U)

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
        if problem.type == "RP":
            ro_l = problem.ro_l
            ro_r = problem.ro_r
            p_l = problem.p_l
            p_r = problem.p_r
            for i in range(i_min, i_max):
                for j in range(j_min, j_max):
                    for k in range(k_min, k_max):
                        if problem.dir=='x':
                            u_l = problem.u_l
                            u_r = problem.u_r
                            v_l = 0.
                            w_l = 0.
                            e_l = eos.gete(ro_l, p_l)
                            E_l = e_l + u_l*u_l/2. + v_l*v_l/2. + w_l*w_l/2.
                            v_r = 0.
                            w_r = 0.
                            e_r = eos.gete(ro_r, p_r)
                            E_r = e_r + u_r*u_r/2. + v_r*v_r/2. + w_r*w_r/2.
                            if self.x_mesh[i] < problem.q_0 and math.fabs(self.x_mesh[i]-problem.q_0)>self.dx/100.:
                                self.U[i][j][k] = [ro_l, ro_l*u_l, ro_l*v_l, ro_l*w_l, ro_l*E_l]
                            else:
                                self.U[i][j][k] = [ro_r, ro_r*u_r, ro_r*v_r, ro_r*w_r, ro_r*E_r]
                        elif problem.dir == 'y':
                            u_l = 0.
                            v_l = problem.u_l
                            w_l = 0.
                            e_l = eos.gete(ro_l, p_l)
                            E_l = e_l + u_l * u_l / 2. + v_l * v_l / 2. + w_l * w_l / 2.
                            u_r = 0.
                            v_r = problem.u_r
                            w_r = 0.
                            e_r = eos.gete(ro_r, p_r)
                            E_r = e_r + u_r * u_r / 2. + v_r * v_r / 2. + w_r * w_r / 2.
                            if self.y_mesh[j] < problem.q_0 and math.fabs(self.y_mesh[j] - problem.q_0) > self.dy / 100.:
                                self.U[i][j][k] = [ro_l, ro_l * u_l, ro_l * v_l, ro_l * w_l, ro_l * E_l]
                            else:
                                self.U[i][j][k] = [ro_r, ro_r * u_r, ro_r * v_r, ro_r * w_r, ro_r * E_r]
                        elif problem.dir == 'z':
                            u_l = 0.
                            v_l = 0.
                            w_l = problem.u_l
                            e_l = eos.gete(ro_l, p_l)
                            E_l = e_l + u_l * u_l / 2. + v_l * v_l / 2. + w_l * w_l / 2.
                            u_r = 0.
                            v_r = 0.
                            w_r = problem.u_r
                            e_r = eos.gete(ro_r, p_r)
                            E_r = e_r + u_r * u_r / 2. + v_r * v_r / 2. + w_r * w_r / 2.
                            if self.z_mesh[k] < problem.q_0 and math.fabs(self.z_mesh[k] - problem.q_0) > self.dz / 100.:
                                self.U[i][j][k] = [ro_l, ro_l * u_l, ro_l * v_l, ro_l * w_l, ro_l * E_l]
                            else:
                                self.U[i][j][k] = [ro_r, ro_r * u_r, ro_r * v_r, ro_r * w_r, ro_r * E_r]
                        else:
                            print("Error: CField.set_ic(): Sorry, only x-direction case can be considered. Bye!")
                            exit(-1)
        elif problem.type == "RTI":
            U = self.U
            ro_down = problem.ro_down
            ro_up = problem.ro_up
            u = 0.
            v = 0.
            w = 0.
            p_0 = problem.p_0
            g = problem.g
            q_0 = problem.q_0
            p = 0.
            for i in range(i_min, i_max):
                for j in range(j_min, j_max):
                    for k in range(k_min, k_max):
                        x = .5*self.dx + self.x_mesh[i]
                        y = .5*self.dy + self.y_mesh[j]
                        z = .5*self.dz + self.z_mesh[k]
                        if problem.dir == 'x':
                            q = x
                        elif problem.dir == 'y':
                            q = y
                        else:
                            q = z
                        if q < q_0:
                            ro = ro_down
                        else:
                            ro = ro_up
                        p = p_0 + ro*g*(q - q_0)
                        e = eos.gete(ro, p)
                        E = e + .5*(0.*0. + 0.*0. + 0.*0.)
                        self.U[i][j][k] = [ro, ro*u, ro*v, ro*w, ro*E]
            # Apply initial disturbance
            # Uncomment the variant you prefer
            # Yalinewich 2D disturbance
            PI = 3.14159
            w_0 = 0.0025
            for i in range(i_min, i_max):
                for j in range(j_min, j_max):
                    for k in range(k_min, k_max):
                        x = self.dx * (.5 + self.x_mesh[i])
                        y = self.dy * (.5 + self.y_mesh[j])
                        z = self.dz * (.5 + self.z_mesh[k])
                        if problem.dir == 'x':
                            self.U[i][j][k][3] = 0.
                            self.U[i][j][k][1] = self.U[i][j][k][0]*w_0* \
                                                    (1. - math.cos(4.*PI*z)) * (1.-math.cos(4.*PI*x/3.))
                        elif problem.dir == 'y':
                            U[i][j][k][1] = 0.
                            U[i][j][k][2] = U[i][j][k][0]*w_0*(1. - math.cos(4.*PI*x)) * (1.-math.cos(4.*PI*y/3.))
                        elif problem.dir == 'z':
                            self.U[i][j][k][2] = 0.
                            self.U[i][j][k][3] = self.U[i][j][k][0]*w_0* \
                                                    (1. - math.cos(4.*PI*y)) * (1.-math.cos(4.*PI*z/3.))
        else:
            print("Error: CField.set_ic(): unknown problem type! Only 1d-PRs and 2d-RTIs allowed. Bye!")
            exit(-1)
        return

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
                    elif bcs[0] == 'w':
                        for num in [0, 2, 3, 4]:  # 0 -> 3, 1 -> 2, i_min-1 -> i_min, i_min-2 -> i_min+1
                            self.U[i][j][k][num] = self.U[self.i_min + (self.i_min - i - 1)][j][k][num]
                        for num in [1]:
                            self.U[i][j][k][num] = - self.U[self.i_min + (self.i_min - i - 1)][j][k][num]
                    else:
                        print("Errof field.set_ics(): only wall-type and transmissive boundaries supported! Bye!")
        # Right X-b.c.
        for i in range(self.i_max, self.i_max+n_bound):
            for j in range(self.j_min, self.j_max):
                for k in range(self.k_min, self.k_max): 
                    if bcs[1] == 't':
                        self.U[i][j][k] = self.U[self.i_max-1][j][k]
                    elif bcs[1] == 'w':
                        for num in [0, 2, 3, 4]:   # i_max -> i_max-1 , i_max+1-> i_max-2
                            self.U[i][j][k][num] = self.U[self.i_max - (i - self.i_max + 1)][j][k][num]
                        for num in [1]:
                            self.U[i][j][k][num] = - self.U[self.i_max - (i - self.i_max + 1)][j][k][num]
                    else:
                        print("Error field.set_ics(): only wall-type and transmissive boundaries supported! Bye!")
        # Left Y-b.c.
        for i in range(0, self.i_max+n_bound):
            for j in range(0, self.j_min):
                for k in range(self.k_min, self.k_max): 
                    if bcs[2] == 't':
                        self.U[i][j][k] = self.U[i][self.j_min][k]
                    elif bcs[2] == 'w':
                        for num in [0, 1, 3, 4]:
                            self.U[i][j][k][num] = self.U[i][self.j_min + (self.j_min - j - 1)][k][num]
                        for num in [2]:
                            self.U[i][j][k][num] = - self.U[i][self.j_min + (self.j_min - j - 1)][k][num]
                    else:
                        print("Error field.set_ics(): only wall-type and transmissive boundaries supported! Bye!")
        # Right Y-b.c.
        for i in range(0, self.i_max+n_bound):
            for j in range(self.j_max, self.j_max+n_bound):
                for k in range(self.k_min, self.k_max): 
                    if bcs[3] == 't':
                        self.U[i][j][k] = self.U[i][self.j_max-1][k]
                    elif bcs[3] == 'w':
                        for num in [0, 1, 3, 4]:
                            self.U[i][j][k][num] = self.U[i][self.j_max - (j - self.j_max + 1)][k][num]
                        for num in [2]:
                            self.U[i][j][k][num] = -self.U[i][self.j_max - (j - self.j_max + 1)][k][num]
                    else:
                        print("Error field.set_ics(): only wall-type and transmissive boundaries supported! Bye!")
        # Left Z-b.c.
        for i in range(0, self.i_max+n_bound):
            for j in range(0, self.j_max+n_bound):
                for k in range(0, self.k_min): 
                    if bcs[4] == 't':
                        self.U[i][j][k] = self.U[i][j][self.k_min]
                    elif bcs[4] == 'w':
                        for num in [0, 1, 2, 4]:
                            self.U[i][j][k][num] = self.U[i][j][self.k_min + (self.k_min - k - 1)][num]
                        for num in [3]:
                            self.U[i][j][k][num] = - self.U[i][j][self.k_min + (self.k_min - k - 1)][num]
                    else:
                        print("Error field.set_ics(): only wall-type and transmissive boundaries supported! Bye!")
        # Right Z-b.c.
        for i in range(0, self.i_max+n_bound):
            for j in range(0, self.j_max+n_bound):
                for k in range(self.k_max, self.k_max+n_bound):
                    if bcs[5] == 't':
                        self.U[i][j][k] = self.U[i][j][self.k_max-1]
                    elif bcs[5] == 'w':
                        for num in [0, 1, 2, 4]:
                            self.U[i][j][k][num] = self.U[i][j][self.k_max - (k - self.k_max + 1)][num]
                        for num in [3]:
                            self.U[i][j][k][num] = - self.U[i][j][self.k_max - (k - self.k_max + 1)][num]
                    else:
                        print("Error field.set_ics(): only wall-type and transmissive boundaries supported! Bye!")

                             
        
                            
                            
                            