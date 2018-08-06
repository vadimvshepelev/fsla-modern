# 'main.py' test program

import numpy as np

import config as cfg
import problem as problem_module
import solver.exact as exc_module
import eos.ideal as eos_module
import field as field_module
import output as output_module
import app as app_module

greeting_str = "================================\n"+ \
               "FSLA eulerian 3D hydrocode v.0.1\n"+ \
               "================================"
print(greeting_str)
eos = eos_module.EOSIdeal(GAMMA=1.4)
# Uncomment for 1d Toro test series
# Don't forget to change 'config' module options

# x:
# problem = [problem_module.CProblem(eos, *problem_module.toro_test_1_x),
#           problem_module.CProblem(eos, *problem_module.toro_test_2_x),
#           problem_module.CProblem(eos, *problem_module.toro_test_3_x),
#           problem_module.CProblem(eos, *problem_module.toro_test_4_x),
#           problem_module.CProblem(eos, *problem_module.toro_test_5_x)]

# y:
#problem = [problem_module.CProblem(eos, *problem_module.toro_test_1_y),
#           problem_module.CProblem(eos, *problem_module.toro_test_2_y),
#           problem_module.CProblem(eos, *problem_module.toro_test_3_y),
#           problem_module.CProblem(eos, *problem_module.toro_test_4_y),
#           problem_module.CProblem(eos, *problem_module.toro_test_5_y)]

# z:
#problem = [problem_module.CProblem(*problem_module.toro_test_1_z),
#           problem_module.CProblem(*problem_module.toro_test_2_z),
#           problem_module.CProblem(*problem_module.toro_test_3_z),
#           problem_module.CProblem(*problem_module.toro_test_4_z),
#           problem_module.CProblem(*problem_module.toro_test_5_z)]
#NX = cfg.const['NX']
#NY = cfg.const['NY']
#NZ = cfg.const['NZ']
#field = [field.CField(problem[0], eos, NX, NY, NZ),
#         field.CField(problem[1], eos, NX, NY, NZ),
#         field.CField(problem[2], eos, NX, NY, NZ),
#         field.CField(problem[3], eos, NX, NY, NZ),
#         field.CField(problem[4], eos, NX, NY, NZ)]
#output = [output_module.COutput(problem[0], eos, field[0]),
#          output_module.COutput(problem[1], eos, field[1]),
#          output_module.COutput(problem[2], eos, field[2]),
#          output_module.COutput(problem[3], eos, field[3]),
#          output_module.COutput(problem[4], eos, field[4])]
#solver = exc_module.CExactRiemannSolver(eos)
# for i in range(5):
#    app = app_module.CApp(problem[i], eos, field[i], solver, output[i])
#    app.run()

problem_RTI = problem_module.CProblem()
problem_RTI.init_RTI(*problem_module.yalinevich_test_y)
NX = cfg.const['NX']
NY = cfg.const['NY']
NZ = cfg.const['NZ']
solver = exc_module.CExactRiemannSolver(eos)
field = field_module.CField(problem_RTI, eos, NX, NY, NZ)
tt_list = [float(i) for i in range(16)]
output = output_module.COutput(problem_RTI, eos, field, tt_list)
app = app_module.CApp(problem_RTI, eos, field, solver, output)
app.run()

