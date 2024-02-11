import os
import numpy as np

from fenics import *
from dolfin_adjoint import *

import matplotlib.pyplot as plt
from matplotlib import cm
from stats import save_dict

from problem import SemilinearProblem, BilinearProblem
from solver import Solver
from simulation_data import SimulationData


data = SimulationData()
solver = Solver()
N = data.N


for Problem in [BilinearProblem, SemilinearProblem]:

    name = Problem().__str__()
    outdir = "output/"+name+"/"
    os.makedirs(outdir, exist_ok=True)
    u_init = Constant(0.0)

    print("\n\n------------------")
    print(name)
    print("------------------\n\n")

    stats = {}

    for n in N:

        print("Discretization parameter n = {}".format(n))
        prob = Problem(n=n, alpha=0.0)

        # Solve problem
        sol = solver(prob, u_init)
        solution_final = sol["control_final"].data
        gradient_final = sol["gradient_final"].data

        # Update statistics
        stats[n] = {"control_final": solution_final.vector()[:],
                    "gradient_final": gradient_final.vector()[:]}

        # Update initial value (homotopy method)
        u_init = solution_final

        # Plot solution
        c = plot(solution_final,wireframe=False, cmap=cm.coolwarm)
        plt.colorbar(c, fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.savefig(outdir+"{}_solution_final_n_{}.pdf".format(name,n))
        plt.close()

    filename = "solutions_gradients"
    save_dict(outdir, filename, stats)
    np.savetxt(outdir  + "/" + filename  + "_filename.txt", np.array([outdir]), fmt = "%s")
