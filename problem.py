class CProblem:
    def __init__(self, GAMMA, dir, ro_l, u_l, p_l, ro_r, u_r, p_r, x_0, t_max):
        """Конструктор одномерной задачи о распаде разрыва в 3D"""
        print("CProblem.__init__():", end="")
        print(" Initializing one-dimensional desintegration of the discontinuity problem...", end="")
        self.direction = dir
        self.ro_l = ro_l
        self.u_l = u_l		
        self.p_l = p_l
        self.ro_r = ro_r
        self.u_r = u_r
        self.p_r = p_r
        self.x_0 = x_0
        self.t_max = t_max
        self.GAMMA = eos_ideal["GAMMA"]
        print("done!")

eos_ideal = {
	"GAMMA":1.4,
}	
	
toro_test_1_x = {
	"direction":'x',	
    "ro_l":1., "u_l":.75, "p_l":1., 
    "ro_r":.125, "u_r":0., "p_r":.1,
	"x_0":.3, "t_max":.2 }
