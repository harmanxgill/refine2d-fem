#!/usr/bin/env python3
"""
Demo script for the 2D FEM with Adaptive Mesh Refinement project.
This script provides an easy way to run the project and see the results.
"""

import sys
import os
import subprocess

# Add the project root to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def main():
    """Main demo function."""
    print("🎯 2D Finite Element Method with Adaptive Mesh Refinement")
    print("=" * 60)
    
    # Check if outputs directory exists, if not create it
    os.makedirs("outputs/vtu", exist_ok=True)
    os.makedirs("outputs/images", exist_ok=True)
    os.makedirs("outputs/animations", exist_ok=True)
    
    print("\n📊 Running adaptive refinement demo...")
    try:
        # Run the adaptive refinement demo
        from scripts.run_adaptive_refinement import run_adaptive_demo
        run_adaptive_demo(nx=8, ny=8, cycles=4, refine_frac=0.3)
        
        print("\n✅ Demo completed successfully!")
        print("\n📁 Output files:")
        print("  • VTU files: outputs/vtu/solution_cycle*.vtu")
        print("  • Images: outputs/images/cycle_*.png")
        print("  • Animation: outputs/animations/refinement_animation.gif")
        
        print("\n🎨 Visualization options:")
        print("  • ParaView: paraview outputs/vtu/solution_cycle0.vtu")
        print("  • Python: python scripts/visualize_solution.py")
        
        print("\n📖 Documentation: docs/README.md")
        
    except Exception as e:
        print(f"\n❌ Error running demo: {e}")
        print("\n💡 Make sure you have installed the dependencies:")
        print("   conda install -c conda-forge numpy scipy meshio matplotlib paraview")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 