import numpy as np
from pycapso.models import get_mean_field


def test_get_mean_field_zero_fixed_point():
    n = 1000

    densities = get_mean_field(num_iter=n, psi0=0, phi0=0)
    result = np.zeros(n)

    assert(np.array_equal(densities[:, 0], result) and
           np.array_equal(densities[:, 1], result))
