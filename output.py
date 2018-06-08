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
        for i in range(self.field.i_min, self.field.i_max):
            for j in range(self.field.j_min, self.field.j_max):
                for k in range(self.field.k_min, self.field.k_max):   
                    f.write("%f %f %f %f %f %f %f %f\n" % (self.field.x_mesh[i], self.field.y_mesh[j], self.field.z_mesh[k], 
                            *self.field.U[i][j][k]))
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
        
    

