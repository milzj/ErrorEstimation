import numpy as np
from fenics import *
from dolfin_adjoint import *
import fw4pde
import moola

class Solver(object):


    def __init__(self):

        maxiter = 100
        gtol = 1e-10
        ftol = -np.inf

        self._options = {"maxiter": maxiter, "gtol": gtol, "ftol": ftol}


    @property
    def options(self):
        return self._options


    def __call__(self, problem, u_init):

        options = self.options

        problem_moola, u_moola = problem(u_init)
        lb = problem.lb
        ub = problem.ub
        beta = problem.beta
        U = problem.control_space
        scaled_L1_norm = problem.scaled_L1_norm

        box_constraints = fw4pde.problem.BoxConstraints(U, lb, ub)
        moola_box_lmo = fw4pde.algorithms.MoolaBoxLMO(box_constraints.lb,
                                                        box_constraints.ub, beta)

        stepsize = fw4pde.stepsize.QuasiArmijoGoldstein(alpha=0.5, gamma=0.75)

        solver = fw4pde.algorithms.FrankWolfe(problem_moola,
                                            initial_point=u_moola,
                                            nonsmooth_functional=scaled_L1_norm,
                                            stepsize=stepsize,
                                            lmo=moola_box_lmo,
                                            options=options)

        return solver.solve()
