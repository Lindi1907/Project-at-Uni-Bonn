import numpy as np
from estimagic import minimize


def sphere(params):
    """Spherical criterion function.

    The unique local and global optimum of this function is at
    the zero vector. It is differentiable, convex and extremely
    well behaved in any possible sense.

    Args:
        params (np.ndarray): 1d numpy array of parameters

    Returns:
        float: The value of the sphere function at params.

    """
    return params @ params


def sphere_variation(params):
    """Less well behaved version of the sphere function."""
    x = params.round(2)
    return x @ x


if __name__ == "__main__":

    res = minimize(
        criterion=sphere,
        params=np.arange(5),
        algorithm="scipy_lbfgsb",
    )
    print(res.params)
