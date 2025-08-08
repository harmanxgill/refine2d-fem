
import numpy as np

def mark_top_fraction(eta, frac=0.3):
    """
    Mark the top fraction of elements by error for refinement.
    Args:
        eta: (M,)
        frac: float
    Returns:
        mask: boolean (M,)
    """
    m = len(eta)
    k = max(1, int(frac * m))
    idx = np.argsort(-eta)[:k]
    mask = np.zeros(m, dtype=bool)
    mask[idx] = True
    return mask

def refine_uniform(coords, tris):
    """
    Placeholder refinement: split every triangle into 4 by edge midpoints.
    Args:
        coords: (N,2)
        tris: (M,3)
    Returns:
        new_coords: (N',2)
        new_tris: (M',3)
    NOTE: This creates hanging nodes; kept as a placeholder for visualization.
    Proper NVB conforming refinement is a TODO.
    """
    # Cache midpoint indices
    edge_to_mid = {}
    def get_mid(i, j):
        a, b = min(i,j), max(i,j)
        key = (a,b)
        if key in edge_to_mid:
            return edge_to_mid[key]
        else:
            mid = 0.5*(coords[a] + coords[b])
            idx = coords.shape[0]
            edge_to_mid[key] = idx
            return idx, mid

    coords_list = [c for c in coords]
    new_tris = []

    for t in tris:
        v0,v1,v2 = t
        # midpoints (we actually have to append to coords)
        m01 = tuple(sorted((v0,v1)))
        m12 = tuple(sorted((v1,v2)))
        m20 = tuple(sorted((v2,v0)))

        for (a,b) in [m01,m12,m20]:
            if (a,b) not in edge_to_mid:
                coords_list.append(0.5*(coords[a]+coords[b]))
                edge_to_mid[(a,b)] = len(coords_list)-1

        i01 = edge_to_mid[m01]
        i12 = edge_to_mid[m12]
        i20 = edge_to_mid[m20]

        # 4 children (red refinement)
        new_tris += [
            [v0, i01, i20],
            [i01, v1, i12],
            [i20, i12, v2],
            [i01, i12, i20],
        ]
    new_coords = np.vstack(coords_list)
    return new_coords, np.array(new_tris, dtype=np.int64)
