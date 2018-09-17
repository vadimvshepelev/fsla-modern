# 'app.py' module, CApp class implementation

import sys, time, math
#from clint.textui import progress

class CApp:
    """Implements the whole calculation process"""
    def __init__(self, problem, eos, field, solver, output):
        print("\nStarting '%s' simulation." % problem.name)
        print("Class CApp: Initializing numerical experiment data...", end="")
        self.problem = problem
        self.eos = eos
        self.field = field
        self.solver = solver	
        self.output = output
        self.t = self.problem.t_min
        self.tau = 0.
        self.counter = 0.
        print("done!")

    def run(self):
        # self.output.write_file(self.problem.name + "-0.dat", self.t)
        self.output.manage_output(self.t)
        init_str = "Calculation in progress:"
        print(init_str)	        
        self.counter = 0
        while self.t < self.problem.t_max:
            self.tau = self.calc_time_step(self.field, self.problem)
            if self.counter < 5:
                self.tau *= .2
            if self.t + self.tau > self.problem.t_max:
                self.tau = self.problem.t_max-self.t
            str_progress_bar = self.output.get_progress_bar(self.t + self.tau)
            output_str = "\r" + str_progress_bar + " iter=%d tau=%.4f t=%.4f CFL=%.1f" % \
                        (self.counter, self.tau, self.t+self.tau, self.problem.CFL)
            sys.stderr.write(output_str)
            self.solver.calc_step(self.field, self.problem, self.tau)

            #if self.t > self.problem.t_max/2 and self.t > self.problem.t_max/2 + self.tau:
            #    self.output.write_file_1d_comp(self.solver, self.problem.name + "-1d-debug.dat", self.t)

            #self.output.write_file_1d_comp(self.solver, self.problem.name + "-1d-debug.dat", self.t)

            #if(self.t > .7*self.problem.t_max):
            #    output_str_file = output_str + " writing to file..."
            #    output_str_file_done = output_str + " writing to file...done!"
            #    output_str_file_clr = output_str + "                        "
            #    sys.stderr.write(output_str_file)
            #    self.output.write_file("test.dat", self.t)
            #    sys.stderr.write(output_str_file_done)
            #    sys.stderr.write(output_str_file_clr)
            self.t += self.tau
            self.counter += 1
            # self.output.write_file(self.problem.name + "-int.dat", self.t)
            self.output.manage_output(self.t)
        print()
        print("done!")
        self.output.write_file(self.problem.name + "-N.dat", self.t)
        self.output.write_file_1d_comp(self.solver, self.problem.name + "-1d.dat", self.t)
        return
    
    def calc_time_step(self, field, problem):
        U = field.U
        i_min = field.i_min
        i_max = field.i_max     
        j_min = field.j_min
        j_max = field.j_max
        k_min = field.k_min
        k_max = field.k_max
        dx = field.dx
        dy = field.dy
        dz = field.dz
        CFL = problem.CFL
        ro0 = U[i_min][j_min][k_min][0]
        u0 = U[i_min][j_min][k_min][1]/ro0
        v0 = U[i_min][j_min][k_min][2]/ro0
        w0 = U[i_min][j_min][k_min][3]/ro0
        e0 = U[i_min][j_min][k_min][4]/ro0 - .5*(u0*u0 + v0*v0 + w0*w0)
        c0 = self.eos.getc(ro0, e0)
#       tau = CFL*1./3*(dx+dy+dz)/c0
        tau = min(dx/(math.fabs(u0)+c0), dy/(math.fabs(v0)+c0), dz/(math.fabs(w0)+c0))
        for i in range(i_min, i_max):
            for j in range(j_min, j_max):
                for k in range(k_min, k_max):
                    ro = U[i][j][k][0]
                    u = U[i][j][k][1]/ro
                    v = U[i][j][k][2]/ro
                    w = U[i][j][k][3]/ro
                    e = U[i][j][k][4]/ro - .5*(u*u + v*v + w*w)
                    c = 0.
                    try:
                        c = self.eos.getc(ro, e)
                    except ValueError:
                        error_file_name = self.problem.name + "-error-dump.dat"
                        print("\nError: CApp.calc_time_step(): negative pressure/density in the node i=%d j=%d k=%d"
                              "Mesh functions written to file %s" % (i, j, k, error_file_name))
                        self.output.write_file(error_file_name, self.t)
                        exit(1)
#                   tau_new = CFL*1./3*(dx+dy+dz)/с
                    tau_new = CFL*min(dx/(math.fabs(u)+c), dy/(math.fabs(v)+c), dz/(math.fabs(w)+c))
                    if tau_new < tau:
                        tau = tau_new
        return tau
                    
        
         
        