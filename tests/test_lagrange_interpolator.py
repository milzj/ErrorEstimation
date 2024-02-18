import pytest

from dolfin import *


def test_lagrange_interpolator():

    mesh = UnitSquareMesh(20, 20)
    U = FunctionSpace(mesh, "DG", 0)

    u_expr = Expression("sin(pi*x[0])*(x[1]+1)", degree=0)
    u = project(u_expr, U)


    mesh1 = UnitSquareMesh(40, 40)
    U1 = FunctionSpace(mesh, "DG", 0)
    u1 = Function(U1)
    u1.interpolate(u)

    _u1 = Function(U1)
    LagrangeInterpolator.interpolate(_u1, u)


    assert errornorm(_u1, u1, degree_rise=0) == 0
