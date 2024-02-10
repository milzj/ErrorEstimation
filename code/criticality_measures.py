
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

    def __init__(self,x,g,lb,ub,beta,tau=1.0):

        self._x = x
        self._g = g
        self._lb = lb
        self._ub = ub
        self._beta = beta
        self._tau = tau

        self.canonical_residual
        self.canonical_map

        self.normal_residual
        self.normal_map

    @property
    def canonical_residual(self):
        """Evaluated the canonical residual

        x - prox_{\psi/tau}(x-(1/tau)*g)

        """

        x = self._x
        g = self._g
        lb = self._lb
        ub = self._ub
        beta = self._beta
        tau = self._tau

        prox_v = prox_box_l1(x-(1/tau)*g, lb, ub, beta/tau)
        self._canonical_residual = x - prox_v

    @property
    def canonical_map(self):
        """Computes the 2-norm of the canonical residual."""
        return np.linalg.norm(self._canonical_residual)

    @property
    def normal_residual(self, v=None):
        """Evaluated the normal residual

        tau*(v-prox(v)) + g(v).

        If v is not supplied, the function computes
        v according to v = x - (1/tau)*g and assumes that
        g(x) = g(prox(v)).
        """

        x = self._x
        g = self._g
        lb = self._lb
        ub = self._ub
        beta = self._beta
        tau = self._tau

        if v == None:
            v = x - (1/tau)*g

        prox_v = prox_box_l1(v, lb, ub, beta/tau)

        self._normal_residual = tau*(v-prox_v)+g

    @property
    def normal_map(self):
        """Computes the 2-norm of the normal residual."""
        return np.linalg.norm( self._normal_residual)
