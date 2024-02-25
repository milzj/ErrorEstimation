from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from figure_style import *

import sys
from problem import SemilinearProblem, BilinearProblem, LinearProblem
from experiment import Experiment
from stats import save_dict, load_dict

cmap = cmap_blue_orange()

now = sys.argv[1]

for Problem in [BilinearProblem, SemilinearProblem]:

    name = Problem().__str__()

    data = Experiment(name)
    N = data.N

    print("\n\n------------------")
    print(name)
    print("------------------\n\n")

    outdir = "output/{}/{}/".format(now,name)
    filename = name + "_solutions_gradients_{}".format(now)

    stats = load_dict(outdir, filename)

    for n in N:

        mesh = UnitSquareMesh(n, n)
        U = FunctionSpace(mesh, "DG", 0)
        u = Function(U)
        u_vec = stats[n]["control_final"]
        u.vector().set_local(u_vec)

        c = plot(u,wireframe=False, cmap=cmap)
        plt.colorbar(c, fraction=0.046, pad=0.04)
        plt.gca().set_box_aspect(1)
        plt.tight_layout()

        for fmt in ["pdf", "png"]:
            plt.savefig(outdir+"{}_solution_final_n_{}_{}.{}".format(name,n,now,fmt), bbox_inches='tight')

        plt.close()
