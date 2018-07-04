# 'main.py' test program

import numpy as np

import config as cfg
import problem as problem_module
import solver.exact as exc_module
import eos.ideal as eos_module
import field as field
import output as output_module
import app as app_module

greeting_str = "================================\n"+ \
               "FSLA eulerian 3D hydrocode v.0.1\n"+ \
               "================================"
print(greeting_str)
eos = eos_module.EOSIdeal(GAMMA=1.4)
problem = [problem_module.CProblem(eos, *problem_module.toro_test_1_x),
           problem_module.CProblem(eos, *problem_module.toro_test_2_x),
           problem_module.CProblem(eos, *problem_module.toro_test_3_x),
           problem_module.CProblem(eos, *problem_module.toro_test_4_x),
           problem_module.CProblem(eos, *problem_module.toro_test_5_x)]
field = [field.CField(problem[0], eos, NX=100, NY=10, NZ=10),
         field.CField(problem[1], eos, NX=200, NY=10, NZ=10),
         field.CField(problem[2], eos, NX=100, NY=10, NZ=10),
         field.CField(problem[3], eos, NX=100, NY=10, NZ=10),
         field.CField(problem[4], eos, NX=100, NY=10, NZ=10)]
output = [output_module.COutput(problem[0], eos, field[0]),
          output_module.COutput(problem[1], eos, field[1]),
          output_module.COutput(problem[2], eos, field[2]),
          output_module.COutput(problem[3], eos, field[3]),
          output_module.COutput(problem[4], eos, field[4])]
solver = exc_module.CExactRiemannSolver(eos)
#for i in range(5):
#    app = app_module.CApp(problem[i], eos, field[i], solver, output[i])
#    app.run()

app = app_module.CApp(problem[1], eos, field[1], solver, output[1])
app.run()