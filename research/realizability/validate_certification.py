import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from fn_complexity import exchange_rate  # noqa: E402
import common as C  # noqa: E402

CYCLE = [(6, 3, 3), (7, 2, 1), (6, 5, 1)]
t0 = time.time()
for a, b in ((CYCLE[0], CYCLE[1]), (CYCLE[1], CYCLE[2]), (CYCLE[2], CYCLE[0])):
    sa, sb = C.Sig.of(a), C.Sig.of(b)
    A, d, A_mp, d_mp, err = C.certified_A_d(sa, sb)
    lab = -math.log(exchange_rate(a, b))
    lba = -math.log(exchange_rate(b, a))
    print(f"{a} -> {b}:  A={A:+.12f}  A_mp={A_mp:+.12f}  |diff|={err:.2e}  "
          f"package A={(lab - lba) / 2:+.12f}   d={d:.12f} d_mp={d_mp:.12f}")
print(f"{time.time() - t0:.1f}s")
