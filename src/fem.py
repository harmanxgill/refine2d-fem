import numpy as np

def tri_area(coords, tri):
    """
    Calculate the area of a triangle.
    Args:
        coords: array of node coordinates
        tri: array of 3 node indices defining the triangle
    Returns:
        Area of the triangle using the determinant formula.
    """
    x0, y0 = coords[tri[0]]
    x1, y1 = coords[tri[1]]
    x2, y2 = coords[tri[2]]

    return 0.5 * np.abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0))

def tri_grad_phi(coords, tri):
    """
    Computes the gradients of linear (P1) basis functions on a triangle
    Args:
        coords: array of node coordinates
        tri: array of 3 node indices defining the triangle
    Returns:
        Gradients as (3,2) for nodes (v0,v1,v2).
    """
    x = coords[tri]

    x0, y0 = x[0]
    x1, y1 = x[1]
    x2, y2 = x[2]

    # Compute the Jacobian matrix
    J = np.array([x1-x0, x2-x0], [y1-y0, y2-y0], dtype=float)

    detJ = np.linalg.det(J)

    # Reference gradients of linear basis: grad phi0 = [-1,-1], phi1=[1,0], phi2=[0,1]
    ref_grads = np.array([[-1., -1.],
                          [ 1.,  0.],
                          [ 0.,  1.]])

    invJT = np.linalg.inv(J).T
    grads = ref_grads @ invJT
    return grads, 0.5*abs(detJ)