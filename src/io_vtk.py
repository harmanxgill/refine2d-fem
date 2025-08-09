import meshio
import numpy as np
import os

def write_vtu(coords, tris, point_data=None, cell_data=None, path="outputs/vtu/solution.vtu"):
    """
    Write a VTU file.
    Args:
        coords: (N,2)
        tris: (M,3)
        point_data: dict of (N,1) arrays
        cell_data: dict of (M,1) arrays
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cells = [ ("triangle", tris.astype(np.int32)) ]
    pd = {} if point_data is None else point_data
    cd = {} if cell_data is None else cell_data
    mesh = meshio.Mesh(points=coords.astype(float), cells=cells, point_data=pd, cell_data=cd)
    mesh.write(path)
