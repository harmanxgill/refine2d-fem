#!/usr/bin/env python3
"""
Run adaptive mesh refinement cycles and generate visualizations.
This demonstrates the "wow factor" of AMR in action.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import meshio
from src.solve import solve_adaptive

def visualize_cycle(coords, tris, u, eta, cycle, save_path="outputs/images"):
    """Create visualization for a single refinement cycle."""
    
    # Create triangulation
    tri = Triangulation(coords[:, 0], coords[:, 1], tris)
    
    # Create figure with subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: Mesh
    axes[0].triplot(tri, 'k-', linewidth=0.5, alpha=0.7)
    axes[0].set_title(f'Mesh - Cycle {cycle}\n({len(coords)} nodes, {len(tris)} elements)')
    axes[0].set_aspect('equal')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    
    # Plot 2: Solution
    im1 = axes[1].tripcolor(tri, u, shading='flat', cmap='viridis')
    axes[1].set_title(f'Solution u - Cycle {cycle}')
    axes[1].set_aspect('equal')
    plt.colorbar(im1, ax=axes[1])
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('y')
    
    # Plot 3: Error indicators
    im2 = axes[2].tripcolor(tri, eta, shading='flat', cmap='hot')
    axes[2].set_title(f'Error Indicators η - Cycle {cycle}')
    axes[2].set_aspect('equal')
    plt.colorbar(im2, ax=axes[2])
    axes[2].set_xlabel('x')
    axes[2].set_ylabel('y')
    
    plt.tight_layout()
    plt.savefig(f'{save_path}/cycle_{cycle:02d}.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Cycle {cycle}: {len(coords)} nodes, {len(tris)} elements, "
          f"u range [{u.min():.6f}, {u.max():.6f}], "
          f"η range [{eta.min():.6f}, {eta.max():.6f}]")

def run_adaptive_demo(nx=8, ny=8, cycles=4, refine_frac=0.3):
    """Run adaptive refinement demo with visualization."""
    
    print(f"Starting adaptive refinement demo:")
    print(f"  Initial mesh: {nx}x{ny} = {nx*ny*2} elements")
    print(f"  Refinement cycles: {cycles}")
    print(f"  Refinement fraction: {refine_frac}")
    print("-" * 50)
    
    # Create output directories
    os.makedirs("outputs/vtu", exist_ok=True)
    os.makedirs("outputs/images", exist_ok=True)
    os.makedirs("outputs/animations", exist_ok=True)
    
    # Run the adaptive solver
    coords, tris, u, eta = solve_adaptive(nx=nx, ny=ny, cycles=cycles, refine_frac=refine_frac)
    
    print("-" * 50)
    print("Final results:")
    print(f"  Final mesh: {len(coords)} nodes, {len(tris)} elements")
    print(f"  Solution range: [{u.min():.6f}, {u.max():.6f}]")
    print(f"  Error indicator range: [{eta.min():.6f}, {eta.max():.6f}]")
    
    # Generate visualizations for each cycle
    print("\nGenerating visualizations...")
    for cycle in range(cycles):
        # Read the VTU file for this cycle
        vtu_file = f"outputs/vtu/solution_cycle{cycle}.vtu"
        if os.path.exists(vtu_file):
            mesh = meshio.read(vtu_file)
            points = mesh.points[:, :2]
            triangles = mesh.cells[0].data
            u_cycle = mesh.point_data['u']
            eta_cycle = mesh.cell_data['eta'][0]
            
            visualize_cycle(points, triangles, u_cycle, eta_cycle, cycle)
    
    print(f"\nVisualizations saved to outputs/images/cycle_XX.png")
    print("You can create a GIF using:")
    print("  convert -delay 100 -loop 0 outputs/images/cycle_*.png outputs/animations/refinement_animation.gif")

if __name__ == "__main__":
    # Run the demo
    run_adaptive_demo(nx=8, ny=8, cycles=4, refine_frac=0.3) 