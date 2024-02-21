from fenics import *
from dolfin_adjoint import *
from fw4pde.problem import ScaledL1Norm
from .problem import Problem
import moola


class LinearProblem(Problem):

    def __init__(self, n=16, alpha=0.0, mpi_comm=MPI.comm_world):

        super().__init__(n=n, alpha=alpha, dim=1, mpi_comm=mpi_comm)

    def __call__(self, u_init, iterative_solver=False):

        n = self._n
        alpha = self._alpha

        lb = self.lb
        ub = self.ub
        beta = self.beta

        yd = self.yd

        U = self.control_space
        V = self.state_space
        bc = self.boundary_conditions

        u = Function(U)
        u.interpolate(u_init)

        y = TrialFunction(V)
        v = TestFunction(V)

        a = inner(grad(y), grad(v)) * dx
        L = u*v*dx

        A, b  = assemble_system(a, L, bc)

        Y = Function(V)

        if iterative_solver == True:
            solver = KrylovSolver(A, "cg")
        else:
            solver = LUSolver(A)

        solver.solve(Y.vector(), b)

        J = assemble(0.5*inner(Y-yd,Y-yd)*dx) + assemble(0.5*Constant(alpha)*u**2*dx)

        control = Control(u)
        rf = ReducedFunctional(J, control)
        problem_moola = MoolaOptimizationProblem(rf)
        u_moola = moola.DolfinPrimalVector(u)

        return problem_moola, u_moola

    def __str__(self):
        return "LinearProblem"

    @property
    def beta(self):
        return 0.001

    @property
    def lb(self):
        return Constant(-1.0)

    @property
    def ub(self):
        mesh = self.mesh
        return Expression("1+0.1*sin(2*pi*x[0])", degree=0, mpi_comm=mesh.mpi_comm())


    @property
    def yd(self):
        mesh = self.mesh
        return Expression("100.0*x[0]*x[0]", degree = 2, mpi_comm=mesh.mpi_comm())



