
import numpy as np
from scipy.sparse.linalg import cg
from .mesh import unit_square_tri_mesh
from .assemble import assemble_poisson
from .boundary import apply_dirichlet
from .error import zz_error_indicators
from .refine import mark_top_fraction, refine_uniform
from .io_vtk import write_vtu

def manufactured_u(x,y):
    """
    Manufactured solution for testing.
    Args:
        x: (N,)
        y: (N,)
    Returns:
        u: (N,)
    """
    return np.sin(np.pi*x) * np.sin(np.pi*y)

def manufactured_f(x,y):
    """
    Manufactured source term for testing.
    Args:
        x: (N,)
        y: (N,)
    Returns:
        f: (N,)
    """
    return 2*(np.pi**2) * np.sin(np.pi*x) * np.sin(np.pi*y)

def solve_adaptive(nx=8, ny=8, cycles=1, refine_frac=0.3):
    """
    Solve Poisson equation with adaptive mesh refinement.
    Args:
        nx: number of subdivisions in x-direction
        ny: number of subdivisions in y-direction
        cycles: number of refinement cycles
        refine_frac: fraction of elements to refine
    """
    coords, tris, bmask = unit_square_tri_mesh(nx, ny)
    for cycle in range(cycles):
        kappa = lambda x,y: 1.0
        A, b = assemble_poisson(coords, tris, kappa, manufactured_f)

        # Dirichlet g = u_exact on boundary
        g = manufactured_u(coords[:,0], coords[:,1])
        A, b = apply_dirichlet(A, b, bmask, g)

        u, info = cg(A, b, rtol=1e-10, maxiter=200)
        if info != 0:
            print("CG did not fully converge, info=", info)

        # error indicators
        eta = zz_error_indicators(coords, tris, u)

        # write VTK
        write_vtu(coords, tris, point_data={"u": u}, cell_data={"eta":[eta]}, path=f"out/solution_cycle{cycle}.vtu")

        # mark and refine (placeholder uniform refinement)
        marked = mark_top_fraction(eta, frac=refine_frac)
        # TODO: replace with conformity-preserving NVB
        coords, tris = refine_uniform(coords, tris)
        # recompute boundary mask (still unit square)
        tol = 1e-12
        bx = (np.abs(coords[:,0])<tol) | (np.abs(coords[:,0]-1.0)<tol)
        by = (np.abs(coords[:,1])<tol) | (np.abs(coords[:,1]-1.0)<tol)
        bmask = bx | by

    return coords, tris, u, eta
