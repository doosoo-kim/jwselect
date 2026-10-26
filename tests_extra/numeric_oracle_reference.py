#!/usr/bin/env python3
"""Developer-only independent numerical oracle for jwselect 0.5.0.

Run this file to regenerate numeric_oracle_fixture.do, or use --check to
verify that its checked-in contents match. Stata users do not need Python.

The oracle fits the unrestricted coordinates of a full canonical mean by
covariance-weighted minimum distance subject to the requested zero restrictions. Its J
statistic is the full covariance-weighted MD residual quadratic form; it is NOT assigned from
the restriction Wald statistic. No jwselect source or production functions
are imported. Dependencies for fixture generation: Python, NumPy, SciPy.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import numpy as np
from scipy.stats import chi2


@dataclass
class Fit:
    mask: np.ndarray
    df: int
    wald: float
    j: float
    bic: float
    cells: np.ndarray
    cell_cov: np.ndarray
    reports: np.ndarray
    report_cov: np.ndarray
    report_loading: np.ndarray


def gls_fit(z, covariance, precision, targets, p, n, mask):
    """Estimate all free mean coordinates using the covariance-weighted MD normal equations."""
    mask = np.asarray(mask, dtype=int)
    m = len(z)
    active = p + np.flatnonzero(mask)
    free = np.setdiff1d(np.arange(m), active)
    design = np.eye(m)[:, free]
    normal = design.T @ precision @ design
    loading = np.linalg.solve(normal, design.T @ precision)
    fitted = design @ (loading @ z)
    residual = z - fitted
    j = float(residual @ precision @ residual)
    if active.size:
        restricted = z[active]
        wald = float(restricted @ np.linalg.solve(
            covariance[np.ix_(active, active)], restricted))
    else:
        # The unrestricted model fits exactly; remove harmless roundoff.
        assert j < 1e-22
        j = 0.0
        wald = 0.0
    cell_loading = design[:p, :] @ loading
    report_loading = targets @ design @ loading
    cell_cov = cell_loading @ covariance @ cell_loading.T
    report_cov = report_loading @ covariance @ report_loading.T
    np.testing.assert_allclose(j, wald, atol=2e-10, rtol=3e-12)
    np.testing.assert_allclose(
        design.T @ precision @ residual, 0, atol=2e-10, rtol=0)
    return Fit(mask, int(mask.sum()), wald, j, j-mask.sum()*np.log(n),
               cell_loading @ z, (cell_cov+cell_cov.T)/2,
               report_loading @ z, (report_cov+report_cov.T)/2,
               report_loading)


def test_row(fit):
    return np.array([[fit.df, fit.wald, fit.j,
                      chi2.sf(fit.wald, fit.df) if fit.df else np.nan]])


def winner(fits):
    order = sorted(range(len(fits)), key=lambda i: fits[i].bic)
    gap = fits[order[1]].bic-fits[order[0]].bic
    # Keep comparisons comfortably away from objective ties.
    assert gap > 1e-5, (gap, fits[order[0]].mask)
    return fits[order[0]], gap


def make_case(number, rng):
    g = 2 if number % 2 else 3
    p = k = 2*g
    m = p+k
    n = (80, 250, 1000)[(number-1) % 3]
    # Dense cross-cohort and post/pre covariance, with a positive diagonal
    # floor so all canonical matrices remain comfortably positive definite.
    factors = rng.normal(size=(m, m+3))
    covariance = factors @ factors.T/(m+3) + .35*np.eye(m)
    covariance *= (.008, .025, .08)[(number-1) % 3]
    sd = np.sqrt(np.diag(covariance))
    z = rng.normal(.8, .5, size=m)
    if number == 1:
        z[p:] = 0.0                  # exact all-restrictions winner
    elif number == 2:
        z[p:] = sd[p:] * 5*np.sqrt(k*np.log(n))  # exact empty winner
    else:
        z[p:] = sd[p:] * np.sqrt(np.log(n)) * rng.normal(
            0, (.55, .9, 1.25, 1.7)[number % 4], size=k)
    precision = np.linalg.solve(covariance, np.eye(m))
    targets = np.zeros((g+1, m))
    for group in range(g):
        targets[group, 2*group:2*group+2] = .5
    targets[-1] = targets[:-1].mean(axis=0)

    # This lexicographic mask order intentionally differs from the Stata
    # Cartesian candidate order; 03_numeric_oracle.do joins rows by mask.
    masks = np.array(list(product((0, 1), repeat=k)), dtype=int)
    fits = [gls_fit(z, covariance, precision, targets, p, n, mask)
            for mask in masks]
    joint, joint_gap = winner(fits)
    cohort_mask = np.zeros(k, dtype=int)
    local_gaps = []
    for group in range(g):
        local_fits = []
        for bits in product((0, 1), repeat=2):
            mask = np.zeros(k, dtype=int)
            mask[2*group:2*group+2] = bits
            local_fits.append(gls_fit(z, covariance, precision, targets,
                                      p, n, mask))
        local_fit, gap = winner(local_fits)
        cohort_mask |= local_fit.mask
        local_gaps.append(gap)
    cohort = gls_fit(z, covariance, precision, targets, p, n, cohort_mask)
    all_fit = fits[-1]

    # Each raw cohort has [g-3,g-2,g-1,g,g+1]. Its g-1 observation is
    # a deterministic zero normalization, including its covariance row.
    groups = 4+3*np.arange(g)
    cells = np.array([(group, time) for group in groups
                      for time in range(group-3, group+2)])
    postmap = np.array([(group, time) for group in groups
                        for time in (group, group+1)])
    hmap = np.array([(group, time) for group in groups
                     for time in (group-3, group-2)])
    raw_z = np.zeros(5*g)
    raw_cov = np.zeros((5*g, 5*g))
    canon_to_raw = np.array([5*j+t for j in range(g) for t in (3, 4)] +
                            [5*j+t for j in range(g) for t in (0, 1)])
    raw_z[canon_to_raw] = z
    raw_cov[np.ix_(canon_to_raw, canon_to_raw)] = covariance
    transform = np.zeros((m, 5*g))
    for row, raw in enumerate(canon_to_raw):
        transform[row, raw] = 1
        transform[row, 5*(raw//5)+2] = -1
    np.testing.assert_array_equal(transform @ raw_z, z)
    np.testing.assert_array_equal(transform @ raw_cov @ transform.T, covariance)
    assert np.linalg.eigvalsh(covariance).min() > .001
    assert np.linalg.cond(covariance) < 100
    assert np.linalg.eigvalsh(raw_cov).min() > -1e-12
    # A full-covariance effect is present: dropping cross-cohort covariance
    # would change at least the canonical matrix, even if a winner agrees.
    assert np.max(np.abs(covariance[:2, 2:4])) > 1e-5

    matrices = {
        "XO_INFO": [[number, g, p, k, n, joint_gap, min(local_gaps)]],
        "XO_B": raw_z[None, :], "XO_V": raw_cov, "XO_C": cells,
        "XO_Z": z[None, :], "XO_VZ": covariance,
        "XO_HMAP": hmap, "XO_POSTMAP": postmap,
        "XO_TARGETS": targets,
        "XO_JMASK": joint.mask[None, :],
        "XO_JB": joint.reports[None, :], "XO_JV": joint.report_cov,
        "XO_JBC": joint.cells[None, :], "XO_JVC": joint.cell_cov,
        "XO_JLOAD": joint.report_loading,
        "XO_CMASK": cohort.mask[None, :],
        "XO_CB": cohort.reports[None, :], "XO_CV": cohort.report_cov,
        "XO_CBC": cohort.cells[None, :], "XO_CVC": cohort.cell_cov,
        "XO_CLOAD": cohort.report_loading,
        "XO_ALLTEST": test_row(all_fit),
        "XO_JTEST": test_row(joint), "XO_CTEST": test_row(cohort),
        # Column order: K mask bits, df, independent Wald, residual J, WBIC.
        "XO_SCORES": np.array([
            [*fit.mask, fit.df, fit.wald, fit.j, fit.bic] for fit in fits]),
    }
    return matrices


def number_text(value):
    value = float(value)
    return "." if np.isnan(value) else format(value, ".17g")


def matrix_text(name, value):
    value = np.asarray(value, dtype=float)
    assert value.ndim == 2
    nr, nc = value.shape
    assert nr > 0 and 0 < nc <= 15
    # A continuation marker changes physical line layout, not expression
    # complexity. Preallocate and fill one row per command so Stata 17 never
    # has to parse a large joined-matrix literal (the largest is 64 x 10).
    lines = [f"        matrix {name} = J({nr},{nc},.)\n"]
    for row_number, row in enumerate(value, 1):
        literal = ", ".join(number_text(x) for x in row)
        command = f"        matrix {name}[{row_number},1] = ({literal})\n"
        assert len(command) < 500
        lines.append(command)
    return "".join(lines)


def render():
    rng = np.random.default_rng(20261024)
    cases = [make_case(number, rng) for number in range(1, 25)]
    different = sum(not np.array_equal(x["XO_JMASK"], x["XO_CMASK"])
                    for x in cases)
    selected_sizes = [int(np.sum(x["XO_JMASK"])) for x in cases]
    assert different >= 3
    assert 0 in selected_sizes and max(selected_sizes) >= 4
    assert any(0 < n < x["XO_INFO"][0][3]
               for n, x in zip(selected_sizes, cases))
    lines = [
        "*! Independent covariance-weighted MD oracle fixtures; generated by numeric_oracle_reference.py.\n",
        "*! Do not edit numbers by hand. Python is not required to run the Stata tests.\n",
        "*! Matrices are filled row by row to keep Stata expressions short.\n",
        "version 17.0\n",
        "args case_number\n",
    ]
    for i, case in enumerate(cases, 1):
        lines.append('    if "`case_number\'" == "'+str(i)+'" {\n')
        for name, value in case.items():
            lines.append(matrix_text(name, value))
        lines.append("        exit\n    }\n")
    lines.append('    display as error "Unknown numeric oracle case: `case_number\'"\n')
    lines.append("    exit 198\n")
    summary = (f"24 cases; {sum(len(x['XO_SCORES']) for x in cases)} joint candidates; "
               f"{different} joint/cohort mask differences; "
               f"selected restriction counts {sorted(set(selected_sizes))}; "
               f"minimum winner gap "
               f"{min(min(x['XO_INFO'][0][5:]) for x in cases):.6g}")
    return "".join(lines), summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="check the existing fixture without rewriting it")
    args = parser.parse_args()
    output, summary = render()
    destination = Path(__file__).with_name("numeric_oracle_fixture.do")
    if args.check:
        assert destination.read_text(encoding="utf-8") == output, "Fixture mismatch"
        print("PASS fixture consistency: " + summary)
    else:
        destination.write_text(output, encoding="utf-8")
        print("Wrote " + destination.name + ": " + summary)


if __name__ == "__main__":
    main()
