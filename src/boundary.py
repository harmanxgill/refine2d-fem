
import numpy as np
from scipy.sparse import csr_matrix

def apply_dirichlet(A, b, boundary_mask, gvals):
    """
    Modify A,b for Dirichlet BC: u=g on boundary nodes.
    Args:
        A: (N,N) stiffness matrix
        b: (N,) load vector
        boundary_mask: boolean (N,)
        gvals: array (N,) with boundary values (ignored for interior).
    Returns:
        A: (N,N) stiffness matrix with Dirichlet BCs applied    
        b: (N,) load vector with Dirichlet BCs applied
    """
    A = A.tolil()
    N = A.shape[0]
    fixed = np.where(boundary_mask)[0]
    free  = np.where(~boundary_mask)[0]

    # Zero rows of fixed, set diagonal to 1, set b = g
    for i in fixed:
        A.rows[i] = [i]
        A.data[i] = [1.0]
        b[i] = gvals[i]

    # Zero columns of fixed for free rows
    for i in free:
        row_cols = A.rows[i]
        row_data = A.data[i]
        for k, col in enumerate(list(row_cols)):
            if col in set(fixed) and col != i:
                # set A[i, col] = 0
                idx = row_cols.index(col)
                row_cols.pop(idx)
                row_data.pop(idx)

    return A.tocsr(), b
