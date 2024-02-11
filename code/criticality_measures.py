
from prox import prox_box_l1
import numpy as np

class CriticalityMeasures(object):
    """Evaluates several criticality measures

    Computes the criticality measures for the optimization problem

    min_x g^T x + beta*norm(x,1) subject to lb <= x <= ub.

    The proximal operator is computed using a composition formula.

    Parameters:
    -----------
        x : ndarray or float
            input array
        lb, ub : ndarray or float
            lower and upper bounds
        g : ndarray or float
            gradient evaluated at x
        beta : float
            regularization parameter, nonnegative
        tau : float
            criticality measure parameter, positive
    """

    def __init__(self,lb,ub,beta,tau=1.0):

        self._lb = lb
        self._ub = ub
        self._beta = beta
        self._tau = tau

    def prox(self,v):

        lb = self._lb
        ub = self._ub
        beta = self._beta
        tau = self._tau

        return prox_box_l1(v, lb, ub, beta/tau)

    def canonical_residual(self, x, g):
        """Evaluated the canonical residual

        x - prox_{\psi/tau}(x-(1/tau)*g(x))

        """

        lb = self._lb
        ub = self._ub
        beta = self._beta
        tau = self._tau

        prox_v = self.prox(x-(1/tau)*g)
        self._canonical_residual = x - prox_v

    def canonical_map(self, x, g):
        """Computes the 2-norm of the canonical residual."""
        self.canonical_residual(x,g)
        return np.linalg.norm(self._canonical_residual)

    def normal_residual(self, v, g):
        """Evaluated the normal residual

        tau*(v-prox(v)) + g(prox(v)).

        If v is not supplied, the function computes
        v according to v = x - (1/tau)*g and assumes that
        g(x) = g(prox(v)).
        """

        lb = self._lb
        ub = self._ub
        beta = self._beta
        tau = self._tau

        prox_v = self.prox(v)
        self._normal_residual = tau*(v-prox_v)+g

    def normal_map(self, v, g):
        """Computes the 2-norm of the normal residual."""
        self.normal_residual(v,g)
        return np.linalg.norm(self._normal_residual)



