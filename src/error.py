
import numpy as np
from .fem import tri_grad_phi, tri_area

def element_grad_u(coords, tris, u):
    """
    Compute constant grad u_h per triangle.
    Args:
        coords: (N,2)
        tris: (M,3)
        u: (N,)
    Returns:
        grads: (M,2)
    """
    grads = np.zeros((tris.shape[0], 2))
    for e, tri in enumerate(tris):
        g, area2 = tri_grad_phi(coords, tri)
        uh = u[tri]
        # grad u_h = sum_i u_i grad phi_i
        grads[e,:] = uh @ g
    return grads

def zz_error_indicators(coords, tris, u):
    """
    Simple ZZ indicator: L2 norm of (grad u_h - recovered nodal grad) over element.
    Here we approximate by: |T| * ||avg(nodal recovered) - grad(u_h)||^2
    Args:
        coords: (N,2)
        tris: (M,3)
        u: (N,)
    Returns:
        eta: (M,)
    """
    m = tris.shape[0]
    grads_T = element_grad_u(coords, tris, u)            # (m,2)

    # recover nodal gradient by area-weighted average
    n = coords.shape[0]
    sum_g = np.zeros((n,2))
    sum_a = np.zeros(n)

    areas = np.array([tri_area(coords, tri) for tri in tris])
    for e, tri in enumerate(tris):
        for v in tri:
            sum_g[v] += grads_T[e] * areas[e]
            sum_a[v] += areas[e]
    nodal_g = np.zeros((n,2))
    nz = sum_a > 0
    nodal_g[nz] = (sum_g[nz].T / sum_a[nz]).T

    # element indicator: compare element grad to nodal average at its vertices
    eta = np.zeros(m)
    for e, tri in enumerate(tris):
        gbar = nodal_g[tri].mean(axis=0)
        diff = gbar - grads_T[e]
        eta[e] = areas[e] * (diff @ diff)
    return eta
