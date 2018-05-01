# 'main.py' test program

import numpy as np

import config as cfg
import problem as problem
import solver.exact as exc
import eos.ideal as eos
import field as field
import app as app

print("FSLA eulerian 3D hydrocode v.0.1")

eos = eos.EOSIdeal(GAMMA=1.4)
problem = problem.CProblem(eos, *problem.toro_test_1_x)	
#print(problem)
		  
field = field.CField(problem, eos, NX=10, NY=10, NZ=10)
#print(field.x_mesh, field.y_mesh, field.z_mesh)

field.write_file("test.dat", 0.)

solver = exc.CExactRiemannSolver(eos)
app = app.CApp(problem, eos, field, solver)
app.run()


#flux = solver.calc_F(0., 0.)


print("Calculation finished! Thank you, bye!")

