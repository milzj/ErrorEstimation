from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from figure_style import *

import sys
from problem import SemilinearProblem, BilinearProblem, LinearProblem
from experiment import Experiment
from stats import save_dict, load_dict


now = sys.argv[1]

for Problem in [LinearProblem]:

    name = Problem().__str__()

    data = Experiment(name)
    N = data.N

    print("\n\n------------------")
    print(name)
    print("------------------\n\n")

    outdir = "output/{}/{}/".format(now,name)
    filename = "solutions_gradients_{}".format(now)

    stats = load_dict(outdir, filename)

    for n in N:

        mesh = UnitIntervalMesh(n)
        U = FunctionSpace(mesh, "DG", 0)
        u = Function(U)
        u_vec = stats[n]["control_final"]
        u.vector().set_local(u_vec)

        c = plot(u)
        plt.tight_layout()
        plt.gca().set_box_aspect(1)

        for fmt in ["pdf", "png"]:
            plt.savefig(outdir+"{}_solution_final_n_{}_{}.{}".format(name,n,now,fmt), bbox_inches='tight')

        plt.close()
