#!/usr/bin/env python3
"""Independent reference checks and deterministic fixtures for jwselect 0.4.0.

Executes Python/NumPy/SciPy, NOT Stata or Mata.  The companion
joint040_smoke.do exercises the public Stata API against these fixed values.
The oracle estimates all unrestricted mean coordinates by covariance-weighted MD and computes
J from its full residual objective; it does not call jwselect code or obtain
J by simply copying the Wald statistic.

Run: python tests/joint040_reference.py
Regenerate: python tests/joint040_reference.py --write-fixture
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.stats import chi2

P = 2
N = 100
MASKS = np.array(list(itertools.product((0, 1), repeat=2)), dtype=int)
TARGETS = np.array([[1., 0., 0., 0.], [0., 1., 0., 0.], [.5, .5, 0., 0.]])
CELLS = np.array([[3, 1], [3, 2], [3, 3], [4, 2], [4, 3], [4, 4]])


def covariance(rho=.9):
    q = .04 * np.array([[1., rho], [rho, 1.]])
    b = np.array([[.3, -.1], [.2, .4]])
    residual = .05 * np.array([[1., .2], [.2, 1.]])
    return np.block([[residual + b @ q @ b.T, b @ q], [q @ b.T, q]])


def constrained_gls(z, v, p, mask, targets=None):
    """Solve a constrained normal-means regression, retaining every free mean."""
    z, v, mask = np.asarray(z), np.asarray(v), np.asarray(mask)
    assert v.shape == (len(z), len(z))
    assert len(mask) == len(z) - p
    active = p + np.flatnonzero(mask)
    free = np.setdiff1d(np.arange(len(z)), active)
    x = np.eye(len(z))[:, free]
    precision = np.linalg.inv(v)
    # Explicit inverse normal equations are independent of the production
    # restriction-submatrix Cholesky shortcut; fixture covariances are well scaled.
    mapping = np.linalg.solve(x.T @ precision @ x, x.T @ precision)
    fitted = x @ mapping @ z
    residual = z - fitted
    objective = float(residual @ precision @ residual)
    if targets is None:
        targets = np.eye(len(z))[:p]
    action = np.atleast_2d(targets)[:, free] @ mapping
    return objective, action @ z, action @ v @ action.T, action


def restriction_wald(z, v, p, mask):
    active = p + np.flatnonzero(mask)
    return 0. if not len(active) else float(z[active] @ np.linalg.solve(v[np.ix_(active, active)], z[active]))


def select_gls(z, v, p, family, targets=None, n=N):
    objectives = np.array([constrained_gls(z, v, p, m, targets)[0] for m in family])
    scores = objectives - np.sum(family, axis=1) * np.log(n)
    chosen = int(np.argmin(scores))
    return chosen, scores, constrained_gls(z, v, p, family[chosen], targets)


def stats(z, v, masks):
    rows = []
    for mask in masks:
        j, estimates, vv, actions = constrained_gls(z, v, P, mask, TARGETS)
        target, a = TARGETS[-1], actions[-1]
        adjustment = float((target-a) @ z)
        gain = float(target @ v @ (target-a))
        fmsc = adjustment**2 + target @ v @ target - 2*gain
        rows.append([sum(mask), j, j-sum(mask)*np.log(N), estimates[-1], np.sqrt(vv[-1,-1]), adjustment, gain, fmsc])
    return np.asarray(rows)


def raw(z, v):
    # z order is post(g=3), post(g=4), pre(g=3), pre(g=4).
    # Include a zero coefficient/covariance for each normalization cell.
    ix = [2, 5, 0, 3]
    b, vr = np.zeros(6), np.zeros((6, 6))
    b[ix] = z
    vr[np.ix_(ix, ix)] = v
    return b[None, :], vr


def table(z, v, mask):
    j = constrained_gls(z, v, P, mask)[0]
    d = sum(mask)
    return np.array([[d, restriction_wald(z, v, P, mask), j, chi2.sf(j, d) if d else np.nan]])


def fixtures():
    z, v = np.array([.8, 1., .5, .5]), covariance()
    b, vr = raw(z, v)
    _, jb, jv, ja = constrained_gls(z, v, P, [1, 1], TARGETS)
    # Truly block-diagonal by cohort after grouping (post_g, pre_g).
    blockv = np.array([[.06, 0, .01, 0], [0, .07, 0, .015], [.01, 0, .04, 0], [0, .015, 0, .04]])
    blockz = np.array([.8, 1., .2, .5])
    blockb, blockrawv = raw(blockz, blockv)
    blockid, _, (_, blockest, blockvcov, _) = select_gls(blockz, blockv, P, MASKS, TARGETS)
    tiez = np.array([.8, 1., .2, .2])
    tieb, tierawv = raw(tiez, blockv)
    onez, onev = z[[0, 2]], v[np.ix_([0, 2], [0, 2])]
    oneb = np.array([[onez[1], 0, onez[0]]])
    onerawv = np.zeros((3, 3))
    onerawv[np.ix_([2, 0], [2, 0])] = onev
    values = {
        'J40_B': b, 'J40_V': vr, 'J40_C': CELLS,
        'J40_SINGULAR_V': raw(z, covariance(1.))[1],
        'J40_Z': z[None, :], 'J40_VZ': v, 'J40_FAMILY': MASKS,
        'J40_JB': jb[None, :], 'J40_JV': jv, 'J40_JA': ja,
        'J40_CB': np.array([[.8, 1., .9]]),
        'J40_JSTATS': stats(z, v, MASKS),
        'J40_ALLTEST': table(z, v, [1, 1]),
        'J40_ONETEST': table(z, v, [1, 0]),
        'J40_PRE': np.array([[3, 1]]),
        'J40_BLOCK_B': blockb, 'J40_BLOCK_V': blockrawv,
        'J40_BLOCK_MASK': MASKS[blockid][None, :],
        'J40_BLOCK_EST': blockest[None, :], 'J40_BLOCK_VCOV': blockvcov,
        'J40_TIE_B': tieb, 'J40_TIE_V': tierawv,
        'J40_ONE_B': oneb, 'J40_ONE_V': onerawv,
        'J40_ONE_C': np.array([[3, 1], [3, 2], [3, 3]]),
    }
    return values


def check_reference():
    z, v = np.array([.8, 1., .5, .5]), covariance()
    assert np.linalg.eigvalsh(v).min() > 0
    local = np.array([0, 0])
    for g in range(2):
        family = np.zeros((2, 2), dtype=int)
        family[1, g] = 1
        winner, _, _ = select_gls(z, v, P, family)
        local += family[winner]
    joint, _, (_, est, _, _) = select_gls(z, v, P, MASKS, TARGETS)
    assert_array_equal(local, [0, 0])
    assert_array_equal(MASKS[joint], [1, 1])
    assert_allclose(est, [.7, .7, .7], atol=1e-13)
    fullwald = restriction_wald(z, v, P, [1, 1])
    localwald = sum(restriction_wald(z, v, P, mask) for mask in ([1, 0], [0, 1]))
    assert abs(fullwald - localwald) > 5
    assert_allclose(fullwald, 125/19, rtol=1e-14)
    fixture = fixtures()
    assert_array_equal(fixture['J40_BLOCK_MASK'], [[1, 0]])
    # Independent MD objective vs Wald over arbitrary dimensions/targets.
    rng = np.random.default_rng(230426)
    comparisons = 0
    for p, k in ((1, 1), (2, 2), (3, 4), (2, 0)):
        for draw in range(20):
            b = rng.normal(size=(p+k, p+k))
            vv = b @ b.T + .5*np.eye(p+k)
            zz = rng.normal(size=p+k)
            for mask in itertools.product((0, 1), repeat=k):
                j, _, _, _ = constrained_gls(zz, vv, p, mask)
                w = restriction_wald(zz, vv, p, mask)
                assert_allclose(j, w, atol=2e-12, rtol=2e-12)
                comparisons += 1
    # Nontrivial block-diagonal factorization of local and joint objectives.
    blockz = np.array([.8, 1., .2, .5])
    blockv = np.array([[.06, 0, .01, 0], [0, .07, 0, .015], [.01, 0, .04, 0], [0, .015, 0, .04]])
    for mask in MASKS:
        j = constrained_gls(blockz, blockv, P, mask)[0]
        separate = sum(constrained_gls(blockz, blockv, P, [mask[0], 0] if g == 0 else [0, mask[1]])[0] for g in range(2))
        assert_allclose(j, separate, atol=1e-12)
    # Subset tests have the same residual objective whether the other free
    # coordinates are retained or the model contains only one cohort.
    for mask in ([0], [1]):
        onej = constrained_gls(z[[0, 2]], v[np.ix_([0, 2], [0, 2])], 1, mask)[0]
        fullj = constrained_gls(z, v, P, [mask[0], 0])[0]
        assert_allclose(onej, fullj, atol=1e-12)
    oneid, _, _ = select_gls(z[[0, 2]], v[np.ix_([0, 2], [0, 2])], 1, np.array([[0], [1]]))
    assert oneid == 0
    # Ties: compare mathematically identical one-restriction Wald scores.
    tieh = np.array([.8, 1., .2, .2])
    assert restriction_wald(tieh, blockv, P, [1, 0]) == restriction_wald(tieh, blockv, P, [0, 1])
    # Empty tests and empty H have no residual restriction and J=0.
    assert_allclose(constrained_gls(z[:P], v[:P, :P], P, [])[0], 0., atol=1e-25)
    print(f'joint040 numerical reference: PASS ({comparisons} independent MD-J/Wald comparisons).')
    print('Correlated fixture: cohort mask [0, 0]; joint mask [1, 1].')
    print(f'Full joint Wald={fullwald:.12g}; sum of marginal Walds={localwald:.12g}.')
    print('These checks do not execute Stata/Mata; run tests/joint040_smoke.do in Stata.')


def fixture_text():
    lines = ['*! Independent MD oracle; generated by tests/joint040_reference.py.',
             '* Do not edit numeric values manually. Matrices leave the active dataset intact.',
             'version 16.1']
    for name, value in fixtures().items():
        a = np.atleast_2d(value)
        def number(x):
            return '.' if np.isnan(x) else format(float(x), '.17g')
        body = ' \\ '.join(', '.join(number(x) for x in row) for row in a)
        lines.append(f'matrix {name} = ({body})')
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-fixture', action='store_true')
    args = parser.parse_args()
    check_reference()
    path = Path(__file__).with_name('joint040_fixture.do')
    expected = fixture_text()
    if args.write_fixture:
        path.write_text(expected, encoding='utf-8')
        print(f'Wrote {path.name}.')
    else:
        # Allow last-bit BLAS/SciPy differences across platforms while checking
        # every committed number, matrix dimension, and fixture name.
        committed = {}
        for line in path.read_text(encoding='utf-8').splitlines():
            if not line.startswith('matrix '):
                continue
            name, body = line[7:].split(' = ', 1)
            assert body.startswith('(') and body.endswith(')')
            rows = body[1:-1].split(chr(92))
            committed[name] = np.array([[np.nan if x.strip() == '.' else float(x)
                                         for x in row.split(',')] for row in rows])
        reference = fixtures()
        assert committed.keys() == reference.keys(), 'Committed fixture names differ.'
        for name, value in reference.items():
            assert committed[name].shape == np.atleast_2d(value).shape, name
            assert_allclose(committed[name], value, rtol=5e-12, atol=5e-12,
                            equal_nan=True, err_msg=name)
        print('Committed Stata fixture matches the independent generator numerically.')


if __name__ == '__main__':
    main()
