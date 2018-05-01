class CProblem:
    def __init__(self, eos, dir, ro_l, u_l, p_l, ro_r, u_r, p_r, q_0, x_min, x_max, y_min, y_max, z_min, z_max, t_min, t_max, bcs):
        """Конструктор одномерной задачи о распаде разрыва в 3D"""
        print("Class CProblem: ", end="")
        print("Initializing one-dimensional desintegration of the discontinuity problem...", end="")
        self.GAMMA = eos.GAMMA
        self.dir = dir
        self.ro_l = ro_l
        self.u_l = u_l		
        self.p_l = p_l
        self.ro_r = ro_r
        self.u_r = u_r
        self.p_r = p_r
        self.q_0 = q_0
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max        
        self.t_min = t_min
        self.t_max = t_max 
        self.bcs = bcs
        """Boundary conditions transcription: 'w' -- wall, 'p' -- periodic, 't' -- transmissive, 
        order (natural): left X-b.c., right X-b.c., left Y-b.c., right Y-b.c., left Z-b.c., right Z-b.c."""        
        print("done!")
    
    def __str__(self):
        return "<%f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s>" % (self.GAMMA, self.dir, 
                                                       self.ro_l, self.u_l, self.p_l, self.ro_r, self.u_r, self.p_r, self.q_0, 
                                                       self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max, 
                                                       self.t_min, self.t_max, self.bcs)
        
toro_test_1_x = ('x', 1., .75, 1., .125, 0., .1, .3, 0., 1., 0., 1., 0., 1., 0., .2, "tttttt")
