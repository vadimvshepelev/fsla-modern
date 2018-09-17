import config as cfg


class CProblem:
    def __init__(self, eos=None,  name=None, dir=None, ro_l=None, u_l=None, p_l=None, ro_r=None, u_r=None, p_r=None, q_0=None,
                 x_min=None, x_max=None, y_min=None, y_max=None, z_min=None, z_max=None, t_min=None, t_max=None,
                 bcs=None):
        """Конструктор одномерной задачи о распаде разрыва в 3D"""
        print("Class CProblem: ", end="")
        print("Initializing...", end="")
        self.eos = eos
        self.name = name
        self.type = "RP"
        self.dir = dir
        self.ro_l = ro_l
        self.u_l = u_l		
        self.p_l = p_l
        self.ro_r = ro_r
        self.u_r = u_r
        self.p_r = p_r
        self.q_0 = q_0
        self.ro_up = -1.
        self.ro_down = -1.
        self.p_0 = -1.
        self.g = -1.
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max        
        self.t_min = t_min
        self.t_max = t_max 
        self.bcs = bcs
        self.CFL = cfg.const['CFL']

        """Boundary conditions transcription: 'w' -- wall, 'p' -- periodic, 't' -- transmissive, 
        order (natural): left X-b.c., right X-b.c., left Y-b.c., right Y-b.c., left Z-b.c., right Z-b.c."""        
        print("done!")

    def init_RTI(self, name, dir, ro_down, ro_up, p_0, g, q_0,
                 x_min, x_max, y_min, y_max, z_min, z_max, t_min, t_max, bcs):
        """Конструктор задачи о неустойчивости Рэлея-Тейлора в 3D"""
        self.name = name
        self.type = "RTI"
        self.dir = dir
        self.ro_down = ro_down
        self.ro_up = ro_up
        self.p_0 = p_0
        self.g = g
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
        self.CFL = cfg.const['CFL']
        self.type = "RTI"

    def __str__(self):
        return "<%f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s>" % (self.GAMMA, self.name, self.dir,
                                                       self.ro_l, self.u_l, self.p_l, self.ro_r, self.u_r, self.p_r, self.q_0, 
                                                       self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max, 
                                                       self.t_min, self.t_max, self.bcs)


# Pre-installed test problems, suitable for code verification

# 1D Toro tests, various Riemann problems
# x-direction
toro_test_1_x = ('toro-1-x', 'x',      1.,       .75,      1.,    .125,        0.,      .1, .3,
                 0., 1., 0., .1, 0., .1, 0., .2,    "tttttt")
toro_test_2_x = ('toro-2-x', 'x',      1.,       -2.,      .4,      1.,        2.,      .4, .5,
                 0., 1., 0., .05, 0., .05, 0., .15, "tttttt")
toro_test_3_x = ('toro-3-x', 'x',      1.,        0.,   1000.,      1.,        0.,     .01, .5,
                 0., 1., 0., .1, 0., .1, 0., .012,  "tttttt")
toro_test_4_x = ('toro-4-x', 'x', 5.99924,   19.5975, 460.894, 5.99242,  -6.19633, 46.0950, .4,
                 0., 1., 0., .1, 0., .1, 0., .035,  "tttttt")
toro_test_5_x = ('toro-5-x', 'x',      1., -19.59745,   1000.,      1., -19.59745,     .01, .8,
                 0., 1., 0., .1, 0., .1, 0., .012,  "tttttt")
# y-direction
toro_test_1_y = ('toro-1-y', 'y',      1.,       .75,      1.,    .125,        0.,      .1, .3,
                 0., .1, 0., 1., 0., .1, 0., .2,  "tttttt")
toro_test_2_y = ('toro-2-y', 'y',      1.,       -2.,      .4,      1.,        2.,      .4, .5,
                 0., .1, 0., 1., 0., .1, 0., .15,  "tttttt")
toro_test_3_y = ('toro-3-y', 'y',      1.,        0.,   1000.,      1.,        0.,     .01, .5,
                 0., .1, 0., 1., 0., .1, 0., .012, "tttttt")
toro_test_4_y = ('toro-4-y', 'y', 5.99924,   19.5975, 460.894, 5.99242,  -6.19633, 46.0950, .4,
                 0., .1, 0., 1., 0., .1, 0., .035, "tttttt")
toro_test_5_y = ('toro-5-y', 'y',      1., -19.59745,   1000.,      1., -19.59745,     .01, .8,
                 0., 1., 0., 1., 0., .1, 0., .012, "tttttt")
# z-direction
toro_test_1_z = ('toro-1-z', 'z',      1.,       .75,      1.,    .125,        0.,      .1, .3,
                 0., .1, 0., .1, 0., 1., 0., .2,  "tttttt")
toro_test_2_z = ('toro-2-z', 'z',      1.,       -2.,      .4,      1.,        2.,      .4, .5,
                 0., .1, 0., .1, 0., 1., 0., .15,  "tttttt")
toro_test_3_z = ('toro-3-z', 'z',      1.,        0.,   1000.,      1.,        0.,     .01, .5,
                 0., .1, 0., .1, 0., 1., 0., .012, "tttttt")
toro_test_4_z = ('toro-4-z', 'z', 5.99924,   19.5975, 460.894, 5.99242,  -6.19633, 46.0950, .4,
                 0., .1, 0., .1, 0., 1., 0., .035, "tttttt")
toro_test_5_z = ('toro-5-z', 'z',      1., -19.59745,   1000.,      1., -19.59745,     .01, .8,
                 0., 1., 0., .1, 0., 1., 0., .012, "tttttt")

# 2D Rayleigh-Taylor instabilities
yalinevich_test_y = ('yalinevich-RTI-y', 'y', 1., 2., 2.5, -.1, .75, 0., .5, 0., 1.5, 0., .1, 0., 15., 'wwwwww')
