from .criticality_measures import CriticalityMeasures
import fenics

class FEniCSCriticalityMeasures(CriticalityMeasures):

    def __init__(self,function_space,x,g,lb,ub,beta,tau=1.0):

        _x = fenics.project(x, function_space)
        _g = fenics.project(g, function_space)
        _lb = fenics.project(lb, function_space)
        _ub = fenics.project(ub, function_space)

        x_vec = _x.vector().get_local()
        g_vec = _g.vector().get_local()
        lb_vec = _lb.vector().get_local()
        ub_vec = _ub.vector().get_local()

        self.function_space = function_space

        super().__init__(x_vec, g_vec, lb_vec, ub_vec, beta, tau=tau)

    @property
    def canonical_map(self):
        """Computes the L2-norm of the canonical residual."""
        v = fenics.Function(self.function_space)
        v.vector()[:] = self._canonical_residual
        return fenics.norm(v, norm_type = "L2")

    @property
    def normal_map(self):
        """Computes the L2-norm of the normal residual."""
        v = fenics.Function(self.function_space)
        v.vector()[:] = self._normal_residual
        return fenics.norm(v, norm_type = "L2")
