#!/usr/bin/env python3
"""
Simple visualization script for VTU files using matplotlib and meshio.
This provides an alternative to ParaView when OpenGL issues occur.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import meshio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

def visualize_vtu(vtu_file):
    """Visualize a VTU file using matplotlib."""
    
    # Read the VTU file
    mesh = meshio.read(vtu_file)
    
    print(f"Mesh info:")
    print(f"  Number of points: {len(mesh.points)}")
    print(f"  Number of cells: {len(mesh.cells[0].data)}")
    print(f"  Available point data: {list(mesh.point_data.keys())}")
    print(f"  Available cell data: {list(mesh.cell_data.keys())}")
    
    # Extract 2D coordinates (ignore z if present)
    points = mesh.points[:, :2]
    triangles = mesh.cells[0].data  # Assuming first cell type is triangles
    
    # Create triangulation
    tri = Triangulation(points[:, 0], points[:, 1], triangles)
    
    # Create figure with subplots
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: Solution field 'u'
    if 'u' in mesh.point_data:
        u = mesh.point_data['u']
        im1 = axes[0].tripcolor(tri, u, shading='flat', cmap='viridis')
        axes[0].set_title('Solution u')
        axes[0].set_aspect('equal')
        plt.colorbar(im1, ax=axes[0])
        print(f"Solution u range: [{u.min():.6f}, {u.max():.6f}]")
    
    # Plot 2: Error indicators 'eta'
    if 'eta' in mesh.cell_data:
        eta = mesh.cell_data['eta'][0]  # First array in the list
        im2 = axes[1].tripcolor(tri, eta, shading='flat', cmap='hot')
        axes[1].set_title('Error Indicators η')
        axes[1].set_aspect('equal')
        plt.colorbar(im2, ax=axes[1])
        print(f"Error indicators η range: [{eta.min():.6f}, {eta.max():.6f}]")
    
    # Add mesh edges
    for ax in axes:
        ax.triplot(tri, 'k-', linewidth=0.5, alpha=0.3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    
    plt.tight_layout()
    plt.savefig('outputs/images/solution_visualization.png', dpi=150, bbox_inches='tight')
    print("Saved visualization as 'outputs/images/solution_visualization.png'")
    plt.show()

if __name__ == "__main__":
    vtu_file = "outputs/vtu/solution_cycle0.vtu"
    if not os.path.exists(vtu_file):
        print(f"Error: VTU file {vtu_file} not found!")
        print("Please run the Poisson solver first: python examples/run_poisson.py")
        sys.exit(1)
    
    visualize_vtu(vtu_file) 