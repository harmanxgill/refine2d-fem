import numpy as np
from scipy.sparse import coo_matrix

from .fem import tri_grad_phi

def assemble_poisson(coords, tris, kappa_fn, f_fn):
    """
    Assemble global stiffness matrix A and load vector b for Poisson.
    Args:
        coords: (N,2)
        tris: (M,3)
        kappa_fn: (x,y)->scalar
        f_fn: (x,y)->scalar
    Returns:
        A: (N,N) stiffness matrix
    """
    n = coords.shape[0]
    rows, cols, data = [], [], []
    b = np.zeros(n)

    # simple 1-point quadrature at triangle centroid (fine for P1 demo)
    for tri in tris:
        grads, area2 = tri_grad_phi(coords, tri)   # area2 = 2*area
        area = area2
        # centroid
        xc = coords[tri].mean(axis=0)
        kappa = float(kappa_fn(xc[0], xc[1]))
        Ke = kappa * area * (grads @ grads.T)

        # assemble stiffness
        for a in range(3):
            ia = tri[a]
            for b_ in range(3):
                ib = tri[b_]
                rows.append(ia); cols.append(ib); data.append(Ke[a,b_])

        # load (use centroid rule too)
        fe = np.array([f_fn(xc[0], xc[1])/3.0]*3) * area
        for a in range(3):
            b[tri[a]] += fe[a]

    A = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    return A, b