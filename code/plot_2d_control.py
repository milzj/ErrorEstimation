from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

import sys
from problem import SemilinearProblem, BilinearProblem, LinearProblem
from experiment import Experiment
from stats import save_dict, load_dict

plt.rcParams.update({
	"text.usetex": True,
	"text.latex.preamble": r"\usepackage{amsfonts}",
	"font.family": "serif",
	"font.serif": "Computer Modern Roman",
	"font.monospace": "Computer Modern Typewriter",
    "figure.figsize": [5.0, 5.0]})

from matplotlib.colors import LinearSegmentedColormap

def cmap_blue_orange():
	"Color map inspired by cm.coolwarm."
	return LinearSegmentedColormap.from_list(name="cmap_BlueOrange",
           colors =["tab:blue", "lightgrey", "tab:orange"], N=256)

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
    filename = "solutions_gradients_{}".format(now)

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
