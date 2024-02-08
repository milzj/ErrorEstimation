import pytest

import numpy as np
from criticality_measures import CriticalityMeasures

atol = 1e-14

@pytest.mark.parametrize("n", [64, 128, 256])
def test_criticality_measures(n):

    lb = -np.ones(n)
    ub = np.ones(n)
    g = np.random.randn(n)
    x = np.zeros(n)
    beta = 1e-3

    # min g^T x + beta norm(x, L1)

    idx = g > beta
    x[idx] = lb[idx]

    idx = g < -beta
    x[idx] = ub[idx]

    cm = CriticalityMeasures(x, g, lb, ub, beta)

    assert cm.normal_map < atol
    assert cm.canonical_map < atol
