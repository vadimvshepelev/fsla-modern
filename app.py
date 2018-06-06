# 'app.py' module, CApp class implementation

import sys, time
#from clint.textui import progress

class CApp:
    """Implements the whole calculation process"""
    def __init__(self, problem, eos, field, solver, output):
        print("Class CApp: Initializing numerical experiment data...", end="")        
        self.problem = problem
        self.eos = eos
        self.field = field
        self.solver = solver	
        self.output = output
        self.t = self.problem.t_min
        self.tau = 0.        
        print("done!")

    def run(self):
        self.output.write_file("test.dat", self.t)  
        init_str = "Starting computational process..."
        print(init_str)	        
        counter = 0             
        while(True):
            self.tau = self.calc_time_step(self.field, self.problem)
            if self.t + self.tau > self.problem.t_max:
                self.tau = self.problem.t_max-self.t
            #sys.stderr.write("%d: tau=%f t=%f CFL=%f \r" % (counter, self.tau, self.t, self.problem.CFL))
            #str_progress_bar = "[.......   ] 87% "
            str_progress_bar = self.output.get_progress_bar(self.t)
            #output_str = "\r" + str_progress_bar + " iter=%d tau=%f t=%f CFL=%.1f " % (counter, self.tau, self.t, self.problem.CFL)
            output_str = str_progress_bar + " iter=%d tau=%f t=%f CFL=%.1f \n" % (counter, self.tau, self.t, self.problem.CFL)
            #if(self.t < .15*self.problem.t_max):
            #    output_str += " writing data file..."
            #    sys.stderr.write(output_str)
            #    self.output.write_file("test.dat", self.t)        
            #    sys.stderr.write("done!")
            #else:
            #    sys.stderr.write(output_str)
            sys.stderr.write(output_str)
            time.sleep(1.)            
            self.t += self.tau
            if self.t == self.problem.t_max:
                break
            counter += 1        
        print(output_str[:-1], end="")    
        print("...done!")
        self.output.write_file("test.dat", self.t)      
        
        
        
        
        
        
        print(self.output.get_progress_bar(self.problem.t_max))
        
        
        
        
                
    
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
        tau = CFL*1./3*(dx+dy+dz)/c0
        for i in range(i_min, i_max):
            for j in range(j_min, j_max):
                for k in range(k_min, k_max):
                    ro = U[i][j][k][0]
                    u = U[i][j][k][1]/ro
                    v = U[i][j][k][2]/ro
                    w = U[i][j][k][3]/ro
                    e = U[i][j][k][4]/ro - .5*(u*u + v*v + w*w)
                    с = self.eos.getc(ro, e)
                    tau_new = CFL*1./3*(dx+dy+dz)/с
                    if(tau_new < tau): 
                       tau = tau_new
        return tau
                    
        
         
        