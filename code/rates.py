import os
import numpy as np

from fenics import *
from dolfin_adjoint import *

import matplotlib.pyplot as plt
from matplotlib import cm
from stats import save_dict, load_dict

from problem import SemilinearProblem, BilinearProblem, LinearProblem
from experiment import Experiment
from fenics_criticality_measures import FEniCSCriticalityMeasures

from prox import prox_box_l1

from convergence_rates import convergence_rates

set_log_level(30)

now = sys.argv[1]
problem = sys.argv[2]

if problem == "LinearProblem":
    Problems = [LinearProblem]
elif problem == "BilinearProblem":
    Problems = [BilinearProblem]
elif problem == "SemilinearProblem":
    Problems = [SemilinearProblem]

for Problem in Problems:

    name = Problem().__str__()

    data = Experiment(name)
    N = data.N
    Nref = data.Nref
    Alpha = data.Alpha
    alpha = Alpha[0]

    print("\n\n------------------")
    print(name)
    print("N_ref={}\n".format(Nref))
    print("------------------\n\n")

    outdir = "output/{}/{}/".format(now,name)
    filename = name + "_" + "solutions_gradients_{}".format(now)

    normal_maps = []
    canonical_maps = []
    gaps = []
    rgaps = []

    stats = load_dict(outdir, filename)

    # Defining reference problem
    reference_problem = Problem(n=Nref, alpha=alpha,mpi_comm=MPI.comm_world)
    lb_ref = reference_problem.lb
    ub_ref = reference_problem.ub
    beta = reference_problem.beta
    U_ref = reference_problem.control_space
    scaled_L1_norm = reference_problem.scaled_L1_norm
    cm = FEniCSCriticalityMeasures(U_ref, lb_ref, ub_ref, beta)
    print("Creating an instance of the reference problem")
    problem_moola_ref, w_moola_ref = reference_problem(Constant(1.0), iterative_solver=True)
    print("Finished: Creating an instance of the reference problem")


    w_href = Function(reference_problem.control_space)
    _vh = Function(reference_problem.control_space)
    v_href = Function(reference_problem.control_space)
    u_init_ref = Function(reference_problem.control_space)

    for n in N:

        print("Discretization parameter n = {}".format(n))

        # Evaluate canonical map
        problem = Problem(n=n, alpha=alpha, mpi_comm=MPI.comm_self)
        u_init = Function(problem.control_space)
        u_vec = stats[n]["control_final"]
        u_init.vector().set_local(u_vec)

        u_init_ref.interpolate(u_init)

        w_moola_ref.zero()
        w_moola_ref.data.assign(u_init_ref)
        w_moola_ref.bump_version()

        problem_moola_ref.obj(w_moola_ref)
        deriv = problem_moola_ref.obj.derivative(w_moola_ref)
        gradient = deriv.primal()

        canonical_maps.append(cm.canonical_map(w_moola_ref.data, gradient.data))

        # Evaluate (regularized) gap function
        w_href.vector().set_local(cm.proj(u_init_ref.vector().get_local()))
        w_moola_ref.zero()
        w_moola_ref.data.assign(w_href)
        w_moola_ref.bump_version()

        problem_moola_ref.obj(w_moola_ref)
        deriv = problem_moola_ref.obj.derivative(w_moola_ref)
        gradient = deriv.primal()

        gaps.append(cm.gap(w_href, gradient, deriv))
        rgaps.append(cm.rgap(w_href, gradient, deriv))

        # Evaluate normal map
        vh = Function(problem.control_space)
        # Compute vh
        vh_vec = stats[n]["control_final"]-stats[n]["gradient_final"]
        vh.vector().set_local(vh_vec)
        v_href.interpolate(vh)
        v_href_vec = v_href.vector().get_local()
        # Compute w = prox(v)
        w_href.vector().set_local(cm.prox(v_href_vec))

        # Evaluate gradient at w
        w_moola_ref.zero()
        w_moola_ref.data.assign(w_href)
        w_moola_ref.bump_version()

        print("Computing gradient")
        problem_moola_ref.obj(w_moola_ref)
        deriv_ref = problem_moola_ref.obj.derivative(w_moola_ref)
        gradient_ref = deriv_ref.primal()
        print("Finished: Computing gradient")

        normal_maps.append(cm.normal_map(v_href, gradient_ref.data))



    stats = {"n": N, "canonical_map": canonical_maps, "normal_map": normal_maps, "gap": gaps, "rgap": rgaps}
    print(stats)

    convergence_rates(canonical_maps, [1/n for n in N])
    convergence_rates(normal_maps, [1/n for n in N])
    convergence_rates(gaps, [1/n for n in N])
    convergence_rates(rgaps, [1/n for n in N])


    filename = "criticality_measures"
    _filename = name + "_" + filename + "_nref_{}_{}".format(Nref,now)
    save_dict(outdir,_filename, stats)
    np.savetxt(outdir  + filename  + "_filename.txt", np.array([outdir +  _filename ]), fmt = "%s")

