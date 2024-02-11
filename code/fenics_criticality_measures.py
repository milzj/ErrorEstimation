from criticality_measures import CriticalityMeasures
import fenics

class FEniCSCriticalityMeasures(CriticalityMeasures):

    def __init__(self,function_space,lb,ub,beta,tau=1.0):

        _lb = fenics.project(lb, function_space)
        _ub = fenics.project(ub, function_space)

        lb_vec = _lb.vector().get_local()
        ub_vec = _ub.vector().get_local()

        self.function_space = function_space

        super().__init__(lb_vec, ub_vec, beta, tau=tau)

    def canonical_map(self, u, g):
        """Computes the L2-norm of the canonical residual."""

        u_vec = u.vector().get_local()
        g_vec = g.vector().get_local()

        self.canonical_residual(u_vec,g_vec)

        v = fenics.Function(self.function_space)
        v.vector()[:] = self._canonical_residual

        return fenics.norm(v, norm_type = "L2")

    def normal_map(self, v, g):
        """Computes the L2-norm of the normal residual."""
        v_vec = v.vector().get_local()
        g_vec = g.vector().get_local()

        self.normal_residual(v_vec, g_vec)

        w = fenics.Function(self.function_space)
        w.vector()[:] = self._normal_residual

        return fenics.norm(w, norm_type = "L2")
