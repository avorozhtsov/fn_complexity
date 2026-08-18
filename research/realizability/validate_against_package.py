import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from fn_complexity import exchange_rate  # noqa: E402
import common as C  # noqa: E402

random.seed(3)


def rnd():
    n = random.randint(2, 7)
    return tuple(sorted((random.randint(1, 40) for _ in range(n)), reverse=True))


worst_d = worst_A = 0.0
cnt = 0
for _ in range(600):
    a, b = rnd(), rnd()
    if max(a) < 2 or max(b) < 2:
        continue
    sa, sb = C.Sig.of(a), C.Sig.of(b)
    d, A = C.d_and_A(sa, sb)
    lab = -math.log(exchange_rate(a, b))
    lba = -math.log(exchange_rate(b, a))
    worst_d = max(worst_d, abs(d - (lab + lba)))
    worst_A = max(worst_A, abs(A - (lab - lba) / 2))
    cnt += 1
print(f"pairs {cnt}   max |d - package| {worst_d:.3e}   max |A - package| {worst_A:.3e}")
