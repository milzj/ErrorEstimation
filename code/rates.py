import os
import numpy as np

from fenics import *
from dolfin_adjoint import *

import matplotlib.pyplot as plt
from matplotlib import cm
from stats import save_dict, load_dict

from problem import SemilinearProblem, BilinearProblem
from simulation_data import SimulationData
from fenics_criticality_measures import FEniCSCriticalityMeasures

from prox import prox_box_l1

from convergence_rates import convergence_rates

data = SimulationData()
N = data.N
Nref = data.Nref


#for Problem in [BilinearProblem, SemilinearProblem]:
for Problem in [SemilinearProblem]:

    name = Problem().__str__()

    print("\n\n------------------")
    print(name)
    print("N_ref={}\n".format(Nref))
    print("------------------\n\n")

    outdir = "output/{}/".format(name)
    filename = "solutions_gradients"

    normal_maps = []
    canonical_maps = []

    stats = load_dict(outdir, filename)


    reference_problem = Problem(n=Nref, alpha=0.0,mpi_comm=MPI.comm_world)
    lb_ref = reference_problem.lb
    ub_ref = reference_problem.ub
    beta = reference_problem.beta
    U_ref = reference_problem.control_space
    scaled_L1_norm = reference_problem.scaled_L1_norm
    cm = FEniCSCriticalityMeasures(U_ref, lb_ref, ub_ref, beta)

    w_href = Function(reference_problem.control_space)
    _vh = Function(reference_problem.control_space)
    v_href = Function(reference_problem.control_space)

    for n in N:

        print("Discretization parameter n = {}".format(n))

        # Evaluate canonical map
        problem = Problem(n=n, alpha=0.0, mpi_comm=MPI.comm_self)
        u_init = Function(problem.control_space)
        u_vec = stats[n]["control_final"]
        u_init.vector().set_local(u_vec)
        problem_moola, u_moola = reference_problem(u_init)
        problem_moola.obj(u_moola)
        gradient = problem_moola.obj.derivative(u_moola).primal()
        canonical_maps.append(cm.canonical_map(u_moola.data, gradient.data))

        # Evaluate normal map
        vh = Function(problem.control_space)
        # Compute vh
        vh_vec = stats[n]["control_final"]-stats[n]["gradient_final"]
        vh.vector().set_local(vh_vec)
        LagrangeInterpolator.interpolate(_vh, vh)
        v_href.interpolate(_vh)
        v_href_vec = v_href.vector().get_local()
        # Compute w = prox(v)
        w_href.vector().set_local(cm.prox(v_href_vec))
        # Evaluate gradient at w
        problem_moola, w_moola = reference_problem(w_href)
        problem_moola.obj(w_moola)
        gradient = problem_moola.obj.derivative(w_moola).primal()

        normal_maps.append(cm.normal_map(v_href, gradient.data))


    convergence_rates(canonical_maps, [1/n for n in N])
    convergence_rates(normal_maps, [1/n for n in N])

    stats = {"n": N, "canonical_map": canonical_maps, "normal_map": normal_maps}
    print(stats)

    filename = "criticality_measures"
    save_dict(outdir, filename, stats)
    np.savetxt(outdir  + "/" + filename  + "_filename.txt", np.array([outdir]), fmt = "%s")

