import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.solve import solve_adaptive

if __name__ == "__main__":
    # Run 1–2 cycles; refinement is a placeholder (non-conforming) in this starter
    solve_adaptive(nx=12, ny=12, cycles=1, refine_frac=0.3)
    print("Wrote VTK to out/. Open solution_cycle0.vtu in ParaView.")
