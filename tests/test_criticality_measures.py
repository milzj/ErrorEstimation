import pytest

import numpy as np
from criticality_measures import CriticalityMeasures

atol = 1e-14

@pytest.mark.parametrize("n", [64, 128, 256])
def test_criticality_measures(n):

    lb = -np.random.rand(n)
    ub = np.random.rand(n)
    g = np.random.randn(n)
    x = np.zeros(n)
    beta = 1e-3

    # min g^T x + beta norm(x, L1)

    idx = g > beta
    x[idx] = lb[idx]

    idx = g < -beta
    x[idx] = ub[idx]

    cm = CriticalityMeasures(lb, ub, beta)
    v = x-g

    assert cm.normal_map(v,g) < atol
    assert cm.canonical_map(x,g) < atol
    assert cm.rgap(x,g) < atol
