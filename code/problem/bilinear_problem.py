from fenics import *
from dolfin_adjoint import *
from fw4pde.problem import ScaledL1Norm
from .problem import Problem
import moola

#set_log_level(30)

class BilinearProblem(Problem):

    def __init__(self, n=16, alpha=0.0, mpi_comm=MPI.comm_world):

        super().__init__(n=n, alpha=alpha, mpi_comm=mpi_comm)

    def __call__(self, u_init, iterative_solver=False):

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

        y = TrialFunction(V)
        v = TestFunction(V)

        a = (inner(grad(y), grad(v)) - y*u * v) * dx
        L = g*v*dx

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
        return "BilinearProblem"

    @property
    def diffusion_coefficient(self):
        mesh = self.mesh
        return Expression("x[1]*x[1]+0.05", degree = 2, mpi_comm=mesh.mpi_comm())

    @property
    def beta(self):
        return 0.0001

    @property
    def lb(self):
        return Constant(0.0)

    @property
    def ub(self):
        mesh = self.mesh
        return Expression('x[0] <= 0.25 ? 0 : -5+20.0*x[0]', degree=0,  mpi_comm=mesh.mpi_comm())

    @property
    def yd(self):
        mesh = self.mesh
        f = Expression("sin(2.0*pi*x[0])*sin(2*pi*x[1])", degree = 0,  mpi_comm=mesh.mpi_comm())
        return Expression("1.0+f", f=f, degree = 1,  mpi_comm=mesh.mpi_comm())

    @property
    def g(self):
        mesh = self.mesh
        return Expression("10*cos(8*pi*x[0])*cos(8*pi*x[1])", degree = 1, mpi_comm=mesh.mpi_comm())



