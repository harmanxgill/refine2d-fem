import numpy as np

def unit_square_tri_mesh(nx: int, ny:int):
    """
    Generate a uniform triangulation of the unit square [0,1]^2.
    Args:
        nx: number of subdivisions in x-direction
        ny: number of subdivisions in y-direction
    Returns:
        coords: (N,2)
        tris:   (M,3) indices into coords
        boundary_nodes: boolean mask (N,)
    """
    xs = np.linspace(0.0, 1.0, nx+1)
    ys = np.linspace(0.0, 1.0, ny+1)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    coords = np.column_stack([X.ravel(), Y.ravel()])

    # Helper function that converts 2D grid indices (i,j) to 1D vertex indices
    def vid(i, j): return i*(ny+1) + j

    tris = []
    for i in range(nx):
        for j in range(ny):
            v00 = vid(i,j)
            v10 = vid(i+1,j)
            v01 = vid(i,j+1)
            v11 = vid(i+1,j+1)
            # two triangles per cell: (v00, v10, v11) and (v00, v11, v01)
            tris.append([v00, v10, v01])
            tris.append([v10, v11, v01])
    tris = np.array(tris)

    # boundary mask: nodes on edges x=0, x=1, y=0, y=1
    tol = 1e-12
    bx = (np.abs(coords[:,0])<tol) | (np.abs(coords[:,0]-1.0)<tol)
    by = (np.abs(coords[:,1])<tol) | (np.abs(coords[:,1]-1.0)<tol)
    bmask = bx | by
    return coords, tris, bmask