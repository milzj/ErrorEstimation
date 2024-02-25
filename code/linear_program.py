import pytest

import numpy as np
from fenics import *
from dolfin_adjoint import *

import moola
import matplotlib.pyplot as plt

set_log_level(30)

from fw4pde.algorithms import FrankWolfe, MoolaBoxLMO
from fw4pde.problem import ScaledL1Norm, BoxConstraints
from fw4pde.stepsize import DemyanovRubinovOptimalStepSize

# Source https://github.com/dolfin-adjoint/pyadjoint/blob/master/pyadjoint/verification.py
def convergence_rates(E_values, eps_values, show=True):
    from numpy import log
    r = []
    for i in range(1, len(eps_values)):
        r.append(log(E_values[i] / E_values[i - 1])
            / log(eps_values[i] / eps_values[i - 1]))
    if show:
        print("Computed convergence rates: {}".format(r))
    return r

def errornormL1(u, uh, mesh=mesh):
    # TODO: Take difference first
    F = abs(u-uh)*dx(mesh, {'quadrature_degre': 5})
    f = assemble(F)
    return f

def solve_problem(n, n_ref,  u_init=None, maxiter=1000, gtol=1e-15, ftol=-np.inf, discrete_gradient=None):

    set_working_tape(Tape())

    beta = 1e-3
    lb = Expression("-10*x[0]", degree = 0)
    ub = Constant(1.0)

    m = UnitIntervalMesh(n)
    mesh = UnitIntervalMesh(n)

    U = FunctionSpace(mesh, "DG", 0)

    scaled_L1_norm = ScaledL1Norm(U,beta)

    u = Function(U)
    if u_init != None:
        u = project(u_init, U)
#        u_vec = u.vector()[:]
#        lb_vec = project(lb, U).vector()[:]
#        ub_vec = project(ub, U).vector()[:]
#        u_vec = np.clip(u_vec, lb_vec, ub_vec)
#        u.vector()[:] = u_vec

    #J = assemble(Constant(1.0)*u*dx)
    #xi = project(Expression("pow(x[0]-.5,4)*x[0]*(x[0]-1.0)", degree = 0), U); J = assemble(xi*u*dx)
    xi = project(Expression("10*sin(2*pi*x[0])", degree = 0), U); J = assemble(xi*u*dx)

    control = Control(u)
    rf = ReducedFunctional(J, control)

    problem = MoolaOptimizationProblem(rf)
    u_moola = moola.DolfinPrimalVector(u)

    box_constraints = BoxConstraints(U, lb, ub)
    moola_box_lmo = MoolaBoxLMO(box_constraints.lb, box_constraints.ub, beta)

    stepsize =  DemyanovRubinovOptimalStepSize()

    options = {"maxiter": maxiter, "gtol": gtol, "ftol": ftol}

    solver = FrankWolfe(problem, initial_point=u_moola, nonsmooth_functional=scaled_L1_norm, \
                stepsize=stepsize, lmo=moola_box_lmo, options=options)

    sol = solver.solve()

    # compute criticality measure
    solution = sol["control_best"]
    problem.obj(solution)
    gradient = problem.obj.derivative(solution).primal()

    gradient_vec = gradient.data.vector()[:]
    solution_vec = solution.data.vector()[:]

    x_vec = solution_vec - gradient_vec
    w_vec = np.clip(x_vec, -beta, beta)
    lb_vec = project(lb, U).vector()[:]
    ub_vec = project(ub, U).vector()[:]
    w_vec = np.clip(x_vec-w_vec, lb_vec, ub_vec)
    w = Function(U)
    w.vector()[:] = w_vec

    canonical_criticality_measure = errornorm(solution.data, w, degree_rise = 0, mesh=mesh)
    normal_criticality_measure = -1.

    # normal map
    if discrete_gradient != None:
        dg = Function(U)
        dg.interpolate(discrete_gradient)
#        dg.interpolate(gradient.data)

        v_vec = solution_vec - dg.vector()[:]
        w_vec = np.clip(v_vec, -beta, beta)
        prox_v_vec = np.clip(v_vec-w_vec, lb_vec, ub_vec)

        prox_v = Function(U)
        prox_v.vector()[:] = prox_v_vec
        prox_v_moola = moola.DolfinPrimalVector(prox_v)

        problem.obj(prox_v_moola)
        gradient = problem.obj.derivative(prox_v_moola).primal()

        w_vec = v_vec - prox_v_vec
        w = Function(U)
        # take minus
        w.vector()[:] = -w_vec

        normal_criticality_measure = errornorm(gradient.data, w, degree_rise = 0, mesh=mesh)


    solution_error = errornormL1(lb, sol["control_final"].data, mesh)
    solution_error = errornorm(sol["control_final"].data, lb, mesh=mesh)

    return sol["control_best"].data, sol["dual_gap"], canonical_criticality_measure,\
               normal_criticality_measure, solution_error, gradient.data



def test_convergence_rate():
    """Code verification for a one-dimensional initial value problem.

    dual_gap(u_h) should converge with rate h

    distance of u_h to true solution should converge with rate h

    canonical criticality measure should converge with rate h

    normal map-based criticality measure should converge with rate h
    """

    n_ref = 2**18
    ns = [2**n for n in range(2,14)]

    solutions = []
    errors = []
    pg_errors = []
    gradients = []

    gtol = 1e-10
    ftol = -np.inf

    for n in ns:
        print(n)
        print("\n")
        solution, dual_gap, canonical_cm, normal_cm, solution_error, gradient \
            = solve_problem(n, n_ref, u_init=None, maxiter=1000, gtol=gtol, ftol=ftol)
        print(canonical_cm)
        solutions.append(solution)
        errors.append(solution_error)
        gradients.append(gradient)

    dual_gaps = []
    canonical_criticality_measures = []
    normal_criticality_measures = []

    # perform one iteration to get access to dual_gap and criticality measures
    for i in range(np.size(ns)):
        solution, dual_gap, canonical_cm, normal_cm, solution_error, gradient \
            = solve_problem(n_ref, n_ref, u_init=solutions[i], maxiter=0, gtol=gtol, ftol=ftol, \
                discrete_gradient=gradients[i])

        dual_gaps.append(dual_gap)
        canonical_criticality_measures.append(canonical_cm)
        normal_criticality_measures.append(normal_cm)

    # Convergence dual gap
    print("dual_gaps={}".format(dual_gaps))
    rates = convergence_rates(dual_gaps, [1.0/n for n in ns])

    #assert np.isclose(np.median(rates), 1.0, atol=0.2)

    X = np.ones((np.size(ns), 2)); X[:, 1] = np.log([1.0/n for n in ns])
    x, residudals, rank, s = np.linalg.lstsq(X, np.log(dual_gaps), rcond=None)
    rate = x[1]
    constant = np.exp(x[0])

    #assert np.isclose(rate, 1.0, atol=0.2)

    # Convergence of solutions
    print("solution error={}".format(errors))
    rates = convergence_rates(errors, [1.0/n for n in ns])

    #assert np.isclose(np.median(rates), 1.0, atol=0.2)

    X = np.ones((np.size(ns), 2)); X[:, 1] = np.log([1.0/n for n in ns])
    x, residudals, rank, s = np.linalg.lstsq(X, np.log(errors), rcond=None)

    rate = x[1]
    constant = np.exp(x[0])

    #assert np.isclose(rate, 1.0, atol=0.2)

    # Convergence of normal map based criticality measure
    print("normal map={}".format(normal_criticality_measures))
    rates = convergence_rates(normal_criticality_measures, [1.0/n for n in ns])

    #assert np.isclose(np.median(rates), 1.0, atol=0.2)

    X = np.ones((np.size(ns), 2)); X[:, 1] = np.log([1.0/n for n in ns])
    x, residudals, rank, s = np.linalg.lstsq(X, np.log(normal_criticality_measures), rcond=None)
    rate = x[1]
    constant = np.exp(x[0])

    #assert np.isclose(rate, 1.0, atol=0.2)

    # Convergence of canonical criticality measure
    print("canonical map={}".format(canonical_criticality_measures))
    rates = convergence_rates(canonical_criticality_measures, [1.0/n for n in ns])

    #assert np.isclose(np.median(rates), 1.0, atol=0.2)

    X = np.ones((np.size(ns), 2)); X[:, 1] = np.log([1.0/n for n in ns])
    x, residudals, rank, s = np.linalg.lstsq(X, np.log(canonical_criticality_measures), rcond=None)
    rate = x[1]
    constant = np.exp(x[0])

    #assert np.isclose(rate, 1.0, atol=0.2)




if __name__ == "__main__":

    test_convergence_rate()
