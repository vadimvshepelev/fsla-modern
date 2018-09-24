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
NX = cfg.const['NX']
NY = cfg.const['NY']
NZ = cfg.const['NZ']


# Test area
#####################################

import ctypes

from ctypes import *
from ctypes.util import find_library

print(windll.kernel32)
print(cdll.msvcrt)
libc = cdll.msvcrt

hello_ctypes = ctypes.CDLL('./c-test/hello.so').hello
hello_ctypes.restype = ctypes.c_int
i = hello_ctypes()
print("It returned", i)



print("!!! done !!!")
exit(1)





#####################################



# Uncomment for 1d Toro test series
# Don't forget to change 'config' module options
# x:
eos = eos_module.EOSIdeal(GAMMA=1.4)
problem = [problem_module.CProblem(eos, *problem_module.toro_test_1_x),
           problem_module.CProblem(eos, *problem_module.toro_test_2_x),
           problem_module.CProblem(eos, *problem_module.toro_test_3_x),
           problem_module.CProblem(eos, *problem_module.toro_test_4_x),
           problem_module.CProblem(eos, *problem_module.toro_test_5_x)]
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

for i in range(5):
    field = field_module.CField(problem[i], eos, NX, NY, NZ)
    solver = exc_module.CExactRiemannSolver(eos)
    tt_list = [problem[i].t_min, problem[i].t_max]
    output = output_module.COutput(problem[0], eos, field, tt_list)
    app = app_module.CApp(problem[i], eos, field, solver, output)
    app.run()


# Uncomment for 2D-RTI test
# solver = exc_module.CExactRiemannSolver(eos)
# problem_RTI = problem_module.CProblem()
# problem_RTI.init_RTI(*problem_module.yalinevich_test_y)
# field = field_module.CField(problem_RTI, eos, NX, NY, NZ)
# tt_list = [float(i) for i in range(16)]
# output = output_module.COutput(problem_RTI, eos, field, tt_list)
# app = app_module.CApp(problem_RTI, eos, field, solver, output)
# app.run()

