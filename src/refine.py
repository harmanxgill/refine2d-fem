
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

def refine_nvb(coords, tris, marked):
    """
    Newest Vertex Bisection (NVB) refinement for conforming meshes.
    Args:
        coords: (N,2) array of coordinates
        tris: (M,3) array of triangles
        marked: (M,) boolean array indicating which triangles to refine
    Returns:
        new_coords: (N',2) updated coordinates
        new_tris: (M',3) updated triangles
    """
    # Initialize data structures
    n_points = len(coords)
    n_tris = len(tris)
    
    # Copy initial data
    new_coords = coords.copy()
    new_tris = tris.copy()
    
    # Track which triangles have been processed
    processed = np.zeros(n_tris, dtype=bool)
    
    # Edge bisection cache: (v1, v2) -> midpoint_vertex
    edge_bisections = {}
    
    def bisect_edge(v1, v2):
        """Bisect edge between vertices v1 and v2, return midpoint vertex index."""
        nonlocal new_coords
        edge = tuple(sorted([v1, v2]))
        if edge in edge_bisections:
            return edge_bisections[edge]
        
        # Create new midpoint
        midpoint = 0.5 * (new_coords[v1] + new_coords[v2])
        mid_idx = len(new_coords)
        new_coords = np.vstack([new_coords, midpoint])
        edge_bisections[edge] = mid_idx
        return mid_idx
    
    def refine_triangle(tri_idx):
        """Refine a single triangle using NVB."""
        nonlocal new_coords, new_tris, processed
        if tri_idx >= len(processed) or processed[tri_idx]:
            return
        
        tri = new_tris[tri_idx]
        v0, v1, v2 = tri
        
        # Find the longest edge (opposite to newest vertex)
        # In NVB, the newest vertex is typically the last one
        # We'll use the longest edge as the refinement edge
        edges = [
            (np.linalg.norm(new_coords[v1] - new_coords[v2]), v1, v2, v0),
            (np.linalg.norm(new_coords[v2] - new_coords[v0]), v2, v0, v1),
            (np.linalg.norm(new_coords[v0] - new_coords[v1]), v0, v1, v2)
        ]
        edges.sort(reverse=True)  # Sort by length, longest first
        
        # Use the longest edge for bisection
        _, v_long1, v_long2, v_opp = edges[0]
        
        # Bisect the longest edge
        mid_idx = bisect_edge(v_long1, v_long2)
        
        # Create two new triangles
        new_tri1 = [v_long1, mid_idx, v_opp]
        new_tri2 = [mid_idx, v_long2, v_opp]
        
        # Replace the original triangle with the first new one
        new_tris[tri_idx] = new_tri1
        
        # Add the second new triangle
        new_tris = np.vstack([new_tris, new_tri2])
        
        # Extend the processed array for the new triangle
        processed = np.append(processed, False)
        
        # Mark as processed
        processed[tri_idx] = True
        
        # Check if we need to refine the neighbor triangle that shares the bisected edge
        # This is crucial for maintaining conformity
        for other_idx in range(len(new_tris)):
            if other_idx != tri_idx and other_idx < len(processed) and not processed[other_idx]:
                other_tri = new_tris[other_idx]
                if v_long1 in other_tri and v_long2 in other_tri:
                    # This triangle shares the bisected edge, so it must also be refined
                    refine_triangle(other_idx)
                    break
    
    # Refine all marked triangles
    for i in range(n_tris):
        if marked[i]:
            refine_triangle(i)
    
    return new_coords, new_tris

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
