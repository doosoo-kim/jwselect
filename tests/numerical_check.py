#!/usr/bin/env python3
"""Independent mathematical checks for jwselect; this does NOT execute Mata.

Run from any directory: python tests/numerical_check.py
Use --write-fixture to regenerate reference_fixture.json and reference_fixture.do.
Requires numpy and scipy. All matrices here contain sampling covariances.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
V = np.array([[.09, .04, .03, .01],
              [.04, .16, .01, .06],
              [.03, .01, .04, .012],
              [.01, .06, .012, .09]])
Z = np.array([1., 2., -.27, .45])
T = np.eye(4)[:2]
MASKS = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])


def loading(target, mask, covariance=V, p=2):
    """Solve the restricted GLS problem, independently of the Mata shortcut.

    Unrestricted mean coordinates include all posts and unused pre coordinates.
    The selected pre means are fixed to zero. This obtains the target loading
    through GLS normal equations rather than conditional-covariance algebra.
    """
    fixed = p + np.flatnonzero(mask)
    free = np.setdiff1d(np.arange(len(covariance)), fixed)
    x = np.eye(len(covariance))[:, free]
    precision_x = np.linalg.solve(covariance, x)
    beta_map = np.linalg.solve(x.T @ precision_x, precision_x.T)
    return target[free] @ beta_map


def score(target, mask, z=Z, covariance=V, p=2):
    a = loading(target, mask, covariance, p)
    adjustment = (target - a) @ z
    gain = target @ covariance @ (target - a)
    return adjustment**2 + target @ covariance @ target - 2 * gain


def raw_contrasts(cells, baseline_by_group):
    """Build baseline post and optional clean-pre contrasts by cell keys."""
    cells = np.asarray(cells, dtype=int)
    keys = [tuple(row) for row in cells]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate raw cell")
    lookup = {key: i for i, key in enumerate(keys)}
    rows, postmap, premap = [], [], []
    for is_post in (True, False):
        for group in sorted(baseline_by_group):
            base = baseline_by_group[group]
            if (group, base) not in lookup:
                raise ValueError("missing baseline")
            dates = sorted(t for g, t in keys
                           if g == group and (t >= group) == is_post and t != base)
            for date in dates:
                row = np.zeros(len(keys))
                row[lookup[group, date]] = 1
                row[lookup[group, base]] -= 1
                rows.append(row)
                (postmap if is_post else premap).append((group, date))
    return np.array(rows), np.array(postmap), np.array(premap)


def reference_fixture():
    local_masks = np.array([[1, 0], [0, 0]])
    local_a = np.vstack([loading(T[i], local_masks[i]) for i in range(2)])
    global_a = np.vstack([loading(target, [1, 1]) for target in T])
    aggregate_target = T.mean(axis=0)
    a = np.vstack([local_a, local_a.mean(axis=0), loading(aggregate_target, [1, 1])])
    return {
        "description": "Independent GLS oracle; p=2, k=2; V is sampling covariance",
        "N": 500, "z": Z.tolist(), "V": V.tolist(),
        "postmap": [[3, 3], [4, 4]], "hmap": [[3, 1], [4, 1]],
        "groups": [3, 4], "bases": [2, 3], "T": T.tolist(), "w": [.5, .5],
        "family_masks": [[0, 0], [1, 0], [0, 0], [0, 1]],
        "family_owners": [3, 3, 4, 4], "family_ids": [1, 2, 1, 2],
        "focused_both_A": a.tolist(), "focused_both_b": (a @ Z).tolist(),
        "focused_both_V_fixed": (a @ V @ a.T).tolist(),
        "focused_both_selection": [[1, 0], [0, 0], [1, 0], [1, 1]],
        "focused_both_output_map": [[1, 3], [1, 4], [2, 0], [3, 0]],
        "aggregate_fmsc": [score(aggregate_target, mask) for mask in MASKS],
        "wbic_A": np.vstack([global_a, global_a.mean(axis=0)]).tolist(),
        "wbic_b": np.r_[global_a @ Z, global_a.mean(axis=0) @ Z].tolist(),
        "wbic_selection": [[1, 1], [1, 1], [1, 1]],
    }


def check_gls_and_schur():
    """Constrained GLS, conditional means and Schur variance must agree."""
    for mask in MASKS:
        active = 2 + np.flatnonzero(mask)
        a = np.vstack([loading(t, mask) for t in T])
        expected = T.copy()
        schur = V[:2, :2].copy()
        if len(active):
            c = V[:2, active]
            q = V[np.ix_(active, active)]
            expected[:, active] -= np.linalg.solve(q, c.T).T
            schur -= c @ np.linalg.solve(q, c.T)
        assert_allclose(a, expected, atol=2e-15)
        assert_allclose(a @ V @ a.T, schur, atol=2e-15)


def check_shared_covariance():
    """Aggregating standard errors or dropping cross-cohort covariances fails."""
    a = np.vstack([loading(T[0], [1, 0]), loading(T[1], [0, 0])])
    joint = a @ V @ a.T
    w = np.array([.5, .5])
    assert_allclose((w @ a) @ V @ (w @ a), w @ joint @ w)
    assert abs(w @ joint @ w - np.sum(w**2 * np.diag(joint))) > .01


def check_fmsc_unbiased_risk():
    """Analytic expected FMSC equals local MSE, including bias correction."""
    mean = np.array([1., 2., .07, -.11])
    for target, mask in itertools.product(T, MASKS):
        a = loading(target, mask)
        r = target - a
        expected_delta2 = (r @ mean)**2 + r @ V @ r
        expected_score = expected_delta2 + target @ V @ target - 2*target @ V @ r
        true_mse = a @ V @ a + (r @ mean)**2
        assert_allclose(expected_score, true_mse, atol=1e-14)
    # This fixture would choose baseline under the erroneous one-gain score.
    target, mask = T[0], [1, 0]
    a = loading(target, mask)
    delta2, gain = ((target-a) @ Z)**2, target @ V @ (target-a)
    assert delta2 - 2*gain < 0 < delta2 - gain


def check_focused_targets_differ():
    """Cohort selection does not minimize the aggregate-focused criterion."""
    local = [int(np.argmin([score(T[0], [0, 0]), score(T[0], [1, 0])])),
             int(np.argmin([score(T[1], [0, 0]), score(T[1], [0, 1])]))]
    global_mask = MASKS[np.argmin([score(T.mean(axis=0), mask) for mask in MASKS])]
    assert local == [1, 0]
    assert global_mask.tolist() == [1, 1]
    fixture = reference_fixture()
    assert abs(fixture["focused_both_b"][2] - fixture["focused_both_b"][3]) > .1


def check_wbic_joint_projection():
    """Local WBIC choices are assembled before the final JOINT projection."""
    wald = Z[2:]**2 / np.diag(V)[2:]
    assert np.all(wald - np.log(500) < 0)
    joint = np.vstack([loading(t, [1, 1]) for t in T])
    separate = np.vstack([loading(T[0], [1, 0]), loading(T[1], [0, 1])])
    assert np.max(abs(joint @ Z - separate @ Z)) > .01
    assert_allclose(joint @ Z, [1.1963541666666668, 1.7171875])


def raw_fixture():
    cells = np.array([(g, t) for g in (3, 4) for t in range(1, 7)])
    rng = np.random.default_rng(774)
    u = rng.normal(size=(12, 15))
    u *= .15 / np.sqrt(15)
    for i, (g, t) in enumerate(cells):
        if t == g-1:
            u[i] = 0
    covariance = u @ u.T
    z = np.array([.04, 0, .4, .7, .9, 1.1, -.03, .06, 0, .8, 1., 1.4])
    return z, covariance, cells


def check_rebasing_and_order():
    z, covariance, cells = raw_fixture()
    for bases in ({3: 2, 4: 3}, {3: 2, 4: 2}):
        c, post, pre = raw_contrasts(cells, bases)
        perm = np.random.default_rng(11).permutation(len(cells))
        cp, postp, prep = raw_contrasts(cells[perm], bases)
        assert_allclose(c @ z, cp @ z[perm])
        assert_allclose(c @ covariance @ c.T, cp @ covariance[np.ix_(perm, perm)] @ cp.T)
        assert np.array_equal(post, postp) and np.array_equal(pre, prep)
        assert np.linalg.eigvalsh(c @ covariance @ c.T).min() > 0
    common, _, hmap = raw_contrasts(cells, {3: 2, 4: 2})
    # The previously normalized g=4,t=3 zero is now a real pre contrast.
    index = len(common) - len(hmap) + np.flatnonzero(np.all(hmap == [4, 3], axis=1))[0]
    assert_allclose(common[index] @ z, -.06)
    for bad in (np.vstack([cells, cells[0]]), cells[np.any(cells != [3, 2], axis=1)]):
        try:
            raw_contrasts(bad, {3: 2, 4: 3})
        except ValueError:
            pass
        else:
            raise AssertionError("bad raw support was accepted")


def check_focused_envelope():
    """Ellipsoid support bounds cover every candidate and any selected action."""
    rng = np.random.default_rng(99201)
    q = V[2:, 2:]
    radius = np.sqrt(chi2.ppf(.99, 2))
    target = T.mean(axis=0)
    candidates = np.vstack([loading(target, mask) for mask in MASKS])
    corrections = candidates[:, 2:] - target[2:]
    support = radius * np.sqrt(np.einsum("ij,jk,ik->i", corrections, q, corrections))
    for _ in range(2000):
        eps = rng.multivariate_normal(np.zeros(4), V)
        direction = rng.normal(size=2)
        direction *= radius*rng.random()**.5/np.linalg.norm(direction)
        mu = Z[2:] + np.linalg.cholesky(q) @ direction
        actual = candidates @ eps + corrections @ mu
        centered = candidates @ eps + corrections @ Z[2:]
        assert np.all(actual >= centered-support-1e-13)
        assert np.all(actual <= centered+support+1e-13)
        assert actual.min() >= (centered-support).min()-1e-13
        assert actual.max() <= (centered+support).max()+1e-13
    # Component-wise extrema bound all joint choices without enumerating them
    # during inference. Components share the SAME epsilon and nuisance mean.
    components = [[.5*loading(T[0], m) for m in ([0, 0], [1, 0])],
                  [.5*loading(T[1], m) for m in ([0, 0], [0, 1])]]
    for _ in range(500):
        eps = rng.multivariate_normal(np.zeros(4), V)
        direction = rng.normal(size=2)
        direction *= radius*rng.random()**.5/np.linalg.norm(direction)
        mu = Z[2:] + np.linalg.cholesky(q) @ direction
        lows, highs = [], []
        for g, family in enumerate(components):
            a = np.array(family)
            r = a[:, 2:] - .5*T[g, 2:]
            supp = radius*np.sqrt(np.einsum("ij,jk,ik->i", r, q, r))
            center = a @ eps + r @ Z[2:]
            lows.append(np.min(center-supp))
            highs.append(np.max(center+supp))
        for pair in itertools.product(*components):
            a = pair[0]+pair[1]
            actual = a @ eps + (a-target)[2:] @ mu
            assert sum(lows)-1e-13 <= actual <= sum(highs)+1e-13


def check_sampling_scale_and_ties():
    """Selection is unchanged by consistent scaling; first-min ties are stable."""
    target = T.mean(axis=0)
    scores = np.array([score(target, mask) for mask in MASKS])
    c = 13.
    scaled = [score(target, mask, Z*c, V*c*c) for mask in MASKS]
    assert_allclose(scaled, scores*c*c, atol=1e-12)
    assert np.argmin([score(target, [0, 0]), score(target, [0, 0])]) == 0
    # A singular active covariance is not regularized into a different method.
    try:
        np.linalg.cholesky(np.array([[1., 1.], [1., 1.]]))
    except np.linalg.LinAlgError:
        pass
    else:
        raise AssertionError("singular Q was accepted")


def mata_literal(value):
    a = np.atleast_2d(np.asarray(value, dtype=float))
    return "(" + " \\ ".join(", ".join(format(v, ".17g") for v in row) for row in a) + ")"


def write_fixtures():
    f = reference_fixture()
    (HERE / "reference_fixture.json").write_text(json.dumps(f, indent=2)+"\n")
    lines = ["*! Generated by tests/numerical_check.py --write-fixture; do not hand-edit.",
             "*! Requires jwselect Mata functions to be loaded; does not change the dataset.",
             "version 16.1", "capture mata: mata drop jw_test_reference()",
             "mata:", "void jw_test_reference()", "{",
             "struct jw_model scalar JW_M", "struct jw_family scalar JW_F",
             "struct jw_fit scalar JW_R", "real matrix JW_S"]
    for field, value in (("z", np.array(f["z"])[:, None]), ("V", f["V"]),
                         ("postmap", f["postmap"]), ("hmap", f["hmap"]),
                         ("groups", np.array(f["groups"])[:, None]),
                         ("bases", np.array(f["bases"])[:, None]), ("T", f["T"]),
                         ("w", f["w"])):
        lines.append(f"JW_M.{field} = {mata_literal(value)}")
    lines += ["JW_M.p = 2", "JW_M.k = 2", "JW_M.nunit = 500", "JW_M.gap = 1"]
    lines += ["JW_S = jw_infer_sqrt(JW_M.V)",
              "assert(max(abs(JW_S*JW_S'-JW_M.V)) < 1e-12)"]
    for field, value in (("masks", f["family_masks"]),
                         ("owners", np.array(f["family_owners"])[:, None]),
                         ("ids", np.array(f["family_ids"])[:, None])):
        lines.append(f"JW_F.{field} = {mata_literal(value)}")
    lines += ['JW_R = jw_run(JW_M,JW_F,J(0,2,.),"focused","both",4096)']
    checks = [("JW_R.A", f["focused_both_A"]),
              ("JW_R.A*JW_M.z", np.array(f["focused_both_b"])[:, None]),
              ("JW_R.A*JW_M.V*JW_R.A'", f["focused_both_V_fixed"]),
              ("JW_R.selected", f["focused_both_selection"]),
              ("JW_R.outmap", f["focused_both_output_map"])]
    for expression, expected in checks:
        lines.append(f"assert(max(abs({expression} - {mata_literal(expected)})) < 1e-10)")
    lines += ['JW_R = jw_run(JW_M,JW_F,J(0,2,.),"wbic","both",4096)',
              f"assert(max(abs(JW_R.A - {mata_literal(f['wbic_A'])})) < 1e-10)",
              f"assert(max(abs(JW_R.selected - {mata_literal(f['wbic_selection'])})) < 1e-10)",
              "}", "jw_test_reference()", "mata drop jw_test_reference()",
              "end", 'display as result "jwselect Mata reference fixture: PASS"']
    (HERE / "reference_fixture.do").write_text("\n".join(lines)+"\n")
    z, covariance, cells = raw_fixture()
    raw = ["*! Generated raw fixture for smoke.do; sampling covariance has structural zero baselines.",
           f"matrix JW_B = {mata_literal(z)}", f"matrix JW_V = {mata_literal(covariance)}",
           f"matrix JW_C = {mata_literal(cells)}"]
    c, post, pre = raw_contrasts(cells, {3: 2, 4: 3})
    cz, cv = c @ z, c @ covariance @ c.T
    targets = np.zeros((2, len(cz)))
    targets[0, :4], targets[1, 4:7] = .25, 1/3
    for label, mask in (("BASE", [0, 0, 0]), ("MANUAL", [1, 0, 1]),
                        ("WBIC", [1, 1, 1])):
        a = np.vstack([loading(t, mask, cv, 7) for t in targets])
        a = np.vstack([a, a.mean(axis=0)])
        raw += [f"matrix JW_EB_{label} = {mata_literal(a @ cz)}",
                f"matrix JW_EV_{label} = {mata_literal(a @ cv @ a.T)}"]
    (HERE / "raw_fixture.do").write_text("\n".join(raw)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-fixture", action="store_true")
    args = parser.parse_args()
    checks = [check_gls_and_schur, check_shared_covariance, check_fmsc_unbiased_risk,
              check_focused_targets_differ, check_wbic_joint_projection,
              check_rebasing_and_order, check_focused_envelope, check_sampling_scale_and_ties]
    for check in checks:
        check()
        print(f"PASS {check.__name__}")
    if args.write_fixture:
        write_fixtures()
        print("Wrote JSON and Stata reference fixtures")
    print(f"PASS: {len(checks)} independent numerical checks; Stata/Mata NOT executed")


if __name__ == "__main__":
    main()
