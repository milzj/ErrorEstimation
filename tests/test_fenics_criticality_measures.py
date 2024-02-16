import pytest

import numpy as np
from fenics import *
from dolfin_adjoint import *
from fenics_criticality_measures import FEniCSCriticalityMeasures
import moola

atol = 1e-14

@pytest.mark.parametrize("n", [64, 128, 256])
def test_fenics_criticality_measures(n):

    mesh = fenics.UnitSquareMesh(n,n)
    U = fenics.FunctionSpace(mesh, "DG", 0)

    lb = fenics.Constant(-1.0)
    ub = fenics.Constant(1.0)
    g_vec = np.random.randn(U.dim())
    g = fenics.Function(U)
    g.vector()[:] = g_vec


    x = np.zeros(U.dim())
    beta = 1e-3

    # min g^T x + beta norm(x, L1)

    idx = g_vec > beta
    x[idx] = -1.0

    idx = g_vec < -beta
    x[idx] = 1.0

    u = fenics.Function(U)
    u.vector()[:] = x

    cm = FEniCSCriticalityMeasures(U, lb, ub, beta)
    v = fenics.Function(U)
    v.vector()[:] = x - g_vec

    assert cm.normal_map(v,g) < atol
    assert cm.canonical_map(u,g) < atol


@pytest.mark.parametrize("n", [64, 128, 256])
def test_moola_criticality_measures(n):

    mesh = UnitSquareMesh(n,n)
    U = FunctionSpace(mesh, "DG", 0)

    lb = Constant(-1.0)
    ub = Constant(1.0)
    g_vec = np.random.randn(U.dim())
    g = Function(U)
    g.vector()[:] = g_vec

    u = Function(U)
    J = assemble(inner(g, u)*dx)

    control = Control(u)
    rf = ReducedFunctional(J, control)
    problem_moola = MoolaOptimizationProblem(rf)
    u_moola = moola.DolfinPrimalVector(u)
    problem_moola.obj(u_moola)
    deriv = problem_moola.obj.derivative(u_moola)
    grad = deriv.primal()

    x = np.zeros(U.dim())
    beta = 1e-3

    # min g^T x + beta norm(x, L1)

    idx = g_vec > beta
    x[idx] = -1.0

    idx = g_vec < -beta
    x[idx] = 1.0

    u = Function(U)
    u.vector()[:] = x

    cm = FEniCSCriticalityMeasures(U, lb, ub, beta)
    v = Function(U)
    v.vector()[:] = x - g_vec


    assert cm.normal_map(v,g) < atol
    assert cm.canonical_map(u,g) < atol
    assert cm.rgap(u, grad, deriv) < atol
