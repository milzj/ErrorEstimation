from fenics import *
from dolfin_adjoint import *
from fw4pde.problem import ScaledL1Norm

#set_log_level(30)

class Problem(object):

    def __init__(self, n=16, alpha=0.0, mpi_comm=MPI.comm_world):

#        set_working_tape(Tape())

        self._n = n
        self._alpha = alpha

        self._mesh = UnitSquareMesh(mpi_comm,n,n)
        mesh = self.mesh
        self._control_space = FunctionSpace(mesh, "DG", 0)
        self._state_space = FunctionSpace(mesh, "CG", 1)
        V = self.state_space
        self._boundary_conditions = DirichletBC(V, 0.0, "on_boundary")


    def __call__(self):
        raise NotImplementedError("Must be overwritten.")

    @property
    def boundary_conditions(self):
        return self._boundary_conditions

    @property
    def mesh(self):
        return self._mesh

    @property
    def control_space(self):
        return self._control_space

    @property
    def state_space(self):
        return self._state_space

    @property
    def scaled_L1_norm(self):
        U = self.control_space
        beta = self.beta
        return ScaledL1Norm(U,beta)

    @property
    def alpha(self):
        return self._alpha

    @property
    def beta(self):
        raise NotImplementedError("Must be overwritten.")

    @property
    def lb(self):
        raise NotImplementedError("Must be overwritten.")

    @property
    def ub(self):
        raise NotImplementedError("Must be overwritten.")

    @property
    def yd(self):
        raise NotImplementedError("Must be overwritten.")

    @property
    def g(self):
        raise NotImplementedError("Must be overwritten.")

