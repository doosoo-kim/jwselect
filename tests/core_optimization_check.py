#!/usr/bin/env python3
"""Numerical/order checks for the 0.3 core optimizations, NOT execution of Mata.

Requires NumPy and SciPy. Compares multiple-RHS Cholesky solves to independently
computed constrained GLS, including cross-cohort covariance and arbitrary cell
weights. Also checks stable deduplication and preallocated mask ordering.
Run actual Stata smoke tests to verify Mata compilation and end-to-end behavior.
"""
from itertools import product
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.linalg import cho_factor, cho_solve


def gls_loading(v, targets, mask, p):
    """Independent reference: restricted GLS in free mean coordinates."""
    targets = np.atleast_2d(targets)
    active = p + np.flatnonzero(mask)
    free = np.setdiff1d(np.arange(len(v)), active)
    x = np.eye(len(v))[:, free]
    vx = np.linalg.solve(v, x)
    mapping = np.linalg.solve(x.T @ vx, vx.T)
    return targets[:, free] @ mapping


def multiple_rhs(v, z, target, mask, p):
    """New intended numerical path; target loading and Wald solve together."""
    a = target.copy()
    active = p + np.flatnonzero(mask)
    wald = 0.0
    if len(active):
        q = v[np.ix_(active, active)]
        rhs = np.column_stack((v[active] @ target, z[active]))
        solved = cho_solve(cho_factor(q, lower=True), rhs)
        a[active] -= solved[:, 0]
        wald = z[active] @ solved[:, 1]
    adjustment = (target-a) @ z
    gain = target @ v @ (target-a)
    return a, np.array([sum(mask), wald, wald-sum(mask)*np.log(500),
                        a@z, np.sqrt(max(0, a@v@a)), adjustment, gain,
                        adjustment**2 + target@v@target - 2*gain])


def stable_unique_reference(rows, addempty=False):
    result = []
    for row in rows:
        if not any(np.array_equal(row, existing) for existing in result):
            result.append(row.copy())
    k = rows.shape[1]
    if addempty and not any(np.all(row == 0) for row in result):
        result.insert(0, np.zeros(k, dtype=int))
    return np.array(result, dtype=int).reshape(len(result), k)


def stable_unique_sorted(rows, addempty=False):
    n, k = rows.shape
    if not n:
        return np.zeros((int(addempty), k), dtype=int)
    keep = np.zeros(n, dtype=bool)
    if not k:
        keep[0] = True
    else:
        order = sorted(range(n), key=lambda i: tuple(rows[i]))
        first = order[0]
        for previous, current in zip(order[:-1], order[1:]):
            if np.array_equal(rows[current], rows[previous]):
                first = min(first, current)
            else:
                keep[first] = True
                first = current
        keep[first] = True
    out = rows[keep]
    if addempty and not np.any(np.all(out == 0, axis=1)):
        out = np.vstack((np.zeros((1, k), dtype=int), out))
    return out


def check_candidate_order():
    # Ascending canonical dates; bit 1 indexes most recent optional date.
    for k in range(15):
        rows = np.zeros((2**k, k), dtype=int)
        for s in range(1, 2**k):
            for j in range(k):
                rows[s, k-1-j] = (s // 2**j) % 2
        expected = np.array(list(product((0, 1), repeat=k)), dtype=int)
        expected = expected.reshape(2**k, k)
        assert_array_equal(rows, expected)
        assert np.sum(rows[0]) == 0
        assert len(rows) == 2**k
    rng = np.random.default_rng(63001)
    for k in (0, 1, 2, 6):
        for n in (0, 1, 2, 75):
            rows = rng.integers(0, 2, size=(n, k))
            for addempty in (False, True):
                assert_array_equal(stable_unique_sorted(rows, addempty),
                                   stable_unique_reference(rows, addempty))


def check_multiple_rhs():
    rng = np.random.default_rng(63002)
    comparisons = 0
    max_error = 0.0
    for p, k in ((1, 0), (2, 1), (4, 4), (6, 6)):
        m = p+k
        for _ in range(8):
            b = rng.normal(size=(m, m))
            v = b@b.T + np.eye(m)*.3
            z = rng.normal(size=m)
            # Arbitrary post-cell weights; preperiod target entries stay zero.
            targets = np.column_stack((rng.normal(size=(3, p)), np.zeros((3, k))))
            for mask_tuple in product((0, 1), repeat=k):
                mask = np.array(mask_tuple, dtype=int)
                gls = gls_loading(v, targets, mask, p)
                for target, reference in zip(targets, gls):
                    a, stats = multiple_rhs(v, z, target, mask, p)
                    assert_allclose(a, reference, atol=2e-12, rtol=2e-12)
                    active = p+np.flatnonzero(mask)
                    wald = (z[active] @ np.linalg.solve(v[np.ix_(active, active)],
                                                      z[active])) if len(active) else 0
                    r = target-reference
                    expected = [sum(mask), wald, wald-sum(mask)*np.log(500),
                                reference@z, np.sqrt(reference@v@reference),
                                r@z, target@v@r,
                                (r@z)**2 + target@v@target - 2*target@v@r]
                    assert_allclose(stats, expected, atol=2e-10, rtol=2e-11)
                    max_error = max(max_error, np.max(np.abs(stats-expected)))
                    comparisons += 1
                # Joint covariance is retained by the stacked target solve.
                direct = targets.copy()
                if len(active):
                    direct[:, active] -= cho_solve(
                        cho_factor(v[np.ix_(active, active)], lower=True),
                        v[active]@targets.T).T
                assert_allclose(direct@v@direct.T, gls@v@gls.T,
                                atol=2e-10, rtol=2e-11)
    return comparisons, max_error


def check_custom_aggregate():
    rng = np.random.default_rng(63003)
    p, k = 6, 4
    b = rng.normal(size=(p+k, p+k))
    v = b@b.T + np.eye(p+k)
    z = rng.normal(size=p+k)
    w = np.array([.03, .17, .1, .3, .05, .35])
    masks = [np.array([1, 0, 0, 0]), np.array([0, 0, 1, 1])]
    cell_a = np.zeros((p, p+k))
    components = []
    for g, indices in enumerate((range(3), range(3, 6))):
        units = np.eye(p+k)[list(indices)]
        cell_a[list(indices)] = gls_loading(v, units, masks[g], p)
        target = np.zeros(p+k)
        target[list(indices)] = w[list(indices)]
        a, _ = multiple_rhs(v, z, target, masks[g], p)
        components.append(a)
    assert_allclose(w@cell_a, sum(components), atol=2e-12)
    # Unequal within-cohort weights must not become uniformly averaged targets.
    uniform_components = []
    for g, indices in enumerate((range(3), range(3, 6))):
        target = np.zeros(p+k)
        target[list(indices)] = sum(w[list(indices)])/3
        uniform_components.append(multiple_rhs(v, z, target, masks[g], p)[0])
    assert np.max(np.abs(w@cell_a-sum(uniform_components))) > .01


def main():
    check_candidate_order()
    comparisons, max_error = check_multiple_rhs()
    check_custom_aggregate()
    print("PASS: preallocated all-subsets order for K=0..14 and stable deduplication")
    print(f"PASS: {comparisons} multiple-RHS score comparisons with independent GLS")
    print(f"Maximum absolute score difference: {max_error:.3g}")
    print("PASS: full joint covariance and arbitrary cell-weight aggregation")
    print("These are independent Python numerical checks; Mata/Stata was not run.")


if __name__ == '__main__':
    main()
