# 'output.py' module, COutput class implementation

import config as cfg

class COutput:
    """Implements some service functions for file and screen output"""
    def __init__(self, problem, eos, field):
        self.field = field
        self.problem = problem
        self.eos = eos
        self.field = field


    def write_file(self, file_name, t):
        """Dumps mesh function to Tecplot data file for visualization"""
        NX = self.field.i_max-self.field.i_min
        NY = self.field.j_max-self.field.j_min
        NZ = self.field.k_max-self.field.k_min
        # print("Function CField.write_file(): writing U field to file " + file_name + " at time t =", t, "...", end="") 
        f = open(cfg.const["OUTPUT_DIR"]+file_name, 'w')
        f.write('VARIABLES="X","Y","Z","ro","ro*u","ro*v","ro*w","ro*E"\n')
        f.write('TITLE="Conservative variables vector field t = ' + str(t) + '"\n')
        f.write('ZONE T="Numerical", I=%d, J=%d, K=%d, F=POINT\n' % (NX, NY, NZ))
        for k in range(self.field.k_min, self.field.k_max):
            for j in range(self.field.j_min, self.field.j_max):
                for i in range(self.field.i_min, self.field.i_max):
                    f.write("%f %f %f %f %f %f %f %f\n" % (self.field.x_mesh[i], self.field.y_mesh[j], self.field.z_mesh[k],
                            *self.field.U[i][j][k]))
        f.close()        


    def write_file_1d_comp(self, solver, file_name, t):
        """Dumps 1d-slice of mesh function compared with exact solution -- for Riemann problems only"""
        if self.problem.type != "RP":
            print("Error: COutput.write_1d_slice_comp(): non-RP type only RP problems supported!")
            exit(1)
        NX = self.field.i_max - self.field.i_min
        NY = self.field.j_max - self.field.j_min
        NZ = self.field.k_max - self.field.k_min
        NQ = 0
        i = 0
        j = 0
        k = 0
        str_dir = self.problem.dir
        if str_dir == 'x':
            str_vel = 'u'
            NQ = NX
            j = NY // 2
            k = NZ // 2
        elif str_dir == 'y':
            str_vel = 'v'
            NQ = NY
            i = NX // 2
            k = NZ // 2
        elif str_dir == 'z':
            str_vel = 'w'
            NQ = NZ
            i = NX // 2
            j = NY // 2
        else:
            print("Error: COutput.write_1d_slice_comp(): unknown direction symbol, should be 'x', 'y' or 'z'!")
            exit(1)
        f = open(cfg.const["OUTPUT_DIR"] + file_name, 'w')
        f.write('VARIABLES="%s","ro","%s","p","e","ro_ex","%s_ex","p_ex","e_ex"\n' % (str_dir, str_vel, str_vel))
        f.write('TITLE="Primitive variables field t = ' + str(t) + '"\n')
        f.write('ZONE T="Numerical", I=%d F=POINT\n' % NQ)
        if str_dir == 'x':
            for i in range(self.field.i_min, self.field.i_max):
                x = self.field.x_mesh[i]
                x_0 = self.problem.q_0
                ro_l = self.problem.ro_l
                u_l = self.problem.u_l
                p_l = self.problem.p_l
                ro_r = self.problem.ro_r
                u_r = self.problem.u_r
                p_r = self.problem.p_r
                U_ex  = solver.calc_RP_solution(ro_l, u_l, 0., 0., p_l, ro_r, u_r, 0., 0., p_r, x-x_0, t)
                e_ex = self.eos.gete(U_ex.ro, U_ex.p)
                ro = self.field.U[i][j][k][0]
                u = self.field.U[i][j][k][1]/ro
                v = self.field.U[i][j][k][2]/ro
                w = self.field.U[i][j][k][3]/ro
                e = self.field.U[i][j][k][4]/ro - .5*(u*u + v*v + w*w)
                p = self.eos.getp(ro, e)
                f.write("%f %f %f %f %f %f %f %f %f\n" % (x, ro, u, p, e, U_ex.ro, U_ex.u_des, U_ex.p, e_ex))
        elif str_dir == 'y':
            for j in range(self.field.j_min, self.field.j_max):
                y = self.field.y_mesh[j]
                y_0 = self.problem.q_0
                ro_l = self.problem.ro_l
                v_l = self.problem.u_l
                p_l = self.problem.p_l
                ro_r = self.problem.ro_r
                v_r = self.problem.u_r
                p_r = self.problem.p_r
                U_ex  = solver.calc_RP_solution(ro_l, v_l, 0., 0., p_l, ro_r, v_r, 0., 0., p_r, y-y_0, t)
                e_ex = self.eos.gete(U_ex.ro, U_ex.p)
                ro = self.field.U[i][j][k][0]
                u = self.field.U[i][j][k][1]/ro
                v = self.field.U[i][j][k][2]/ro
                w = self.field.U[i][j][k][3]/ro
                e = self.field.U[i][j][k][4]/ro - .5*(u*u + v*v + w*w)
                p = self.eos.getp(ro, e)
                f.write("%f %f %f %f %f %f %f %f %f\n" % (y, ro, v, p, e, U_ex.ro, U_ex.u_des, U_ex.p, e_ex))
        else:
            for k in range(self.field.k_min, self.field.k_max):
                z = self.field.z_mesh[k]
                z_0 = self.problem.q_0
                ro_l = self.problem.ro_l
                w_l = self.problem.u_l
                p_l = self.problem.p_l
                ro_r = self.problem.ro_r
                w_r = self.problem.u_r
                p_r = self.problem.p_r
                U_ex = solver.calc_RP_solution(ro_l, w_l, 0., 0., p_l, ro_r, w_r, 0., 0., p_r, z-z_0, t)
                e_ex = self.eos.gete(U_ex.ro, U_ex.p)
                ro = self.field.U[i][j][k][0]
                u = self.field.U[i][j][k][1]/ro
                v = self.field.U[i][j][k][2]/ro
                w = self.field.U[i][j][k][3]/ro
                e = self.field.U[i][j][k][4]/ro - .5*(u*u + v*v + w*w)
                p = self.eos.getp(ro, e)
                f.write("%f %f %f %f %f %f %f %f %f\n" % (z, ro, w, p, e, U_ex.ro, U_ex.u_des, U_ex.p, e_ex))
        f.close()

    def get_progress_bar(self, t):
        output_string = '['
        progress_rate = (t - self.problem.t_min)/(self.problem.t_max - self.problem.t_min)        
        for counter in range(10+1):
            if int(progress_rate *10) >= counter:
                output_string += '.'
            else:
                output_string += ' '
        output_string += "] %d" % int(progress_rate*100) + '%'
        return output_string
        
    

