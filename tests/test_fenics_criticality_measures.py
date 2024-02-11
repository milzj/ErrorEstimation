import pytest

import numpy as np
import fenics
from fenics_criticality_measures import FEniCSCriticalityMeasures

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
