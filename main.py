# 'main.py' test program

import numpy as np

import config as cfg
import problem as problem
import solver.exact as exc
import eos.ideal as eos
import field as field
import output as output
import app as app

greeting_str = "================================\n"+"FSLA eulerian 3D hydrocode v.0.1\n"+"================================"
print(greeting_str)
eos = eos.EOSIdeal(GAMMA=1.4)
problem = problem.CProblem(eos, *problem.toro_test_1_x)	
#print(problem)
field = field.CField(problem, eos, NX=10, NY=10, NZ=10)
#print(field.x_mesh, field.y_mesh, field.z_mesh)
solver = exc.CExactRiemannSolver(eos)
output = output.COutput(problem, eos, field)
app = app.CApp(problem, eos, field, solver, output)
app.run()


#flux = solver.calc_F(0., 0.)


print("Calculation finished! Thank you, bye!")

