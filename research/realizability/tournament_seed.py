"""Seed measurement for session brief G: which tournaments does the exchange
comparison realise?

Reproduces the numbers quoted in G1.  Run from the repo root.
"""
import itertools, math, random, sys
sys.path.insert(0, "src")
import numpy as np
from fn_complexity import exchange_rate

random.seed(11)


def random_signature():
    n = random.randint(2, 7)
    return tuple(sorted((random.randint(1, 40) for _ in range(n)), reverse=True))


_cache = {}


def L(a, b):
    """L(a,b) = -log C(a->b); the quasi-metric of brief D Part 0."""
    key = (a, b)
    if key not in _cache:
        _cache[key] = -math.log(exchange_rate(a, b))
    return _cache[key]


def flow(signatures):
    """The antisymmetric part A(a,b) = mid_beta(u_a - u_b)."""
    n = len(signatures)
    K = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        v = (L(signatures[i], signatures[j]) - L(signatures[j], signatures[i])) / 2
        K[i, j], K[j, i] = v, -v
    return K


def hodge_split(K):
    """Least-squares gradient part and residual (HodgeRank)."""
    psi = -K.mean(axis=1)
    G = psi[None, :] - psi[:, None]
    nk = np.linalg.norm(K)
    return np.linalg.norm(G) / nk, np.linalg.norm(K - G) / nk


def main():
    pool = list({random_signature() for _ in range(300)})
    print(f"pool: {len(pool)} signatures")

    # G1 -- which 4-vertex tournaments occur?
    seen, used = {}, 0
    for _ in range(1500):
        S = random.sample(pool, 4)
        K = flow(S)
        if np.min(np.abs(K[np.triu_indices(4, 1)])) < 1e-9:
            continue  # a tie: not a tournament
        score = tuple(sorted(int(sum(1 for j in range(4) if K[i, j] > 0)) for i in range(4)))
        seen[score] = seen.get(score, 0) + 1
        used += 1
    print(f"\nG1: 4-vertex score sequences over {used} random 4-subsets")
    print("    (the four isomorphism classes are (0,1,2,3), (1,1,1,3), (0,2,2,2), (1,1,2,2))")
    for k, v in sorted(seen.items()):
        print(f"    {k}  {v}")

    # G3 -- Hodge split and cycle counts at several n
    print("\nG3: Hodge split of the flow A")
    for n in (8, 16, 24):
        S = random.sample(pool, n)
        K = flow(S)
        grad, curl = hodge_split(K)
        cyc = sum(
            1
            for i, j, k in itertools.combinations(range(n), 3)
            if (K[i, j] > 0) == (K[j, k] > 0) == (K[k, i] > 0)
        )
        print(f"    n={n:3d}  |grad|/|A|={grad:.4f}  |curl|/|A|={curl:.4f}  "
              f"3-cycles={cyc} of {math.comb(n,3)}")


if __name__ == "__main__":
    main()
