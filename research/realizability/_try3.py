import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import realize as R  # noqa: E402

# the 3-cycle 0->1->2->0
T = np.zeros((3, 3), dtype=bool)
T[0, 1] = T[1, 2] = T[2, 0] = True
for r in (2, 3, 4):
    t0 = time.time()
    x, m = R.realise(T, r=r, seed=1, maxiter=200, popsize=20)
    ints, p = R.to_integers(R._sigs_from(x, 3, r))
    A, err = R.certified_matrix(ints)
    ok = all((A[i, j] > 0) == T[i, j] for i in range(3) for j in range(3) if i != j)
    cert = min(abs(A[i, j]) for i in range(3) for j in range(3) if i != j)
    print(f"r={r}  margin={m:.3e}  integer margin={cert:.3e}  signs ok={ok}  "
          f"mp err={err:.1e}  {time.time()-t0:.0f}s")
    print("   ", [tuple(s) for s in ints][:3])
