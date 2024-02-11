from fenics import *
from dolfin_adjoint import *
from fw4pde.problem import ScaledL1Norm
from .problem import Problem
import moola

set_log_level(30)


class SemilinearProblem(Problem):

    def __init__(self, n=16, alpha=0.0):

        super().__init__(n=n, alpha=alpha)

    def __call__(self, u_init):

        n = self._n
        alpha = self._alpha

        lb = self.lb
        ub = self.ub
        beta = self.beta

        g = self.g
        yd = self.yd

        U = self.control_space
        V = self.state_space
        bc = self.boundary_conditions

        u = Function(U)
        u.interpolate(u_init)

        y = Function(V)
        v = TestFunction(V)

        F = (inner(grad(y), grad(v)) + y**3 * v - u*v - g*v) * dx
        solve(F == 0, y, bc)

        J = assemble(0.5*inner(y-yd,y-yd)*dx + 0.5*Constant(alpha)*u**2*dx)
        control = Control(u)
        rf = ReducedFunctional(J, control)
        problem_moola = MoolaOptimizationProblem(rf)
        u_moola = moola.DolfinPrimalVector(u)

        return problem_moola, u_moola

    def __str__(self):
        return "SemilinearProblem"

    @property
    def beta(self):
        return 0.005

    @property
    def lb(self):
        return Constant(-10.0)

    @property
    def ub(self):
        return Expression('x[0] <= 0.25 ? 0 : -5.0+20.0*x[0]', degree=0)

    @property
    def yd(self):
        return Expression("2.0*sin(4.0*pi*x[0])*cos(8.0*pi*x[1])*exp(2.0*x[0])", degree = 1)

    @property
    def g(self):
        return Expression("10.0*cos(8.0*pi*x[0])*cos(8.0*pi*x[1])", degree = 1)


