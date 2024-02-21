import os
import numpy as np

from fenics import *
from dolfin_adjoint import *

from stats import save_dict

from problem import SemilinearProblem, BilinearProblem, LinearProblem
from solver import Solver
from experiment import Experiment


set_log_level(30)

solver = Solver()

now = sys.argv[1]

for Problem in [LinearProblem, BilinearProblem, SemilinearProblem]:
#for Problem in [LinearProblem]:

    name = Problem().__str__()
    outdir = "output/"+now+"/"+name+"/"
    os.makedirs(outdir, exist_ok=True)
    u_init = Constant(0.0)

    print("\n\n------------------")
    print(name)
    print("------------------\n\n")

    data = Experiment(name)
    N = data.N
    Nref = data.Nref
    Alpha = data.Alpha
    alpha = Alpha[0]


    stats = {}

    for n in N:

        print("Discretization parameter n = {}".format(n))
        prob = Problem(n=n, alpha=alpha)

        # Solve problem
        sol = solver(prob, u_init)
        solution_final = sol["control_final"].data
        gradient_final = sol["gradient_final"].data

        # Update statistics
        stats[n] = {"control_final": solution_final.vector()[:],
                    "gradient_final": gradient_final.vector()[:]}

        # Update initial value (homotopy method)
        u_init = solution_final

        if "LinearProblem" in name:
            u_init = Constant(0.0)


    filename = "solutions_gradients"
    _filename = name + "_" + filename + "_{}".format(now)
    save_dict(outdir, _filename, stats)
    np.savetxt(outdir + filename  + "_filename.txt", np.array([outdir +"/" + _filename]), fmt = "%s")
