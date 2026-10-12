#!/usr/bin/env python3
"""Independent algebra/scheduling checks for release 0.3.0.

This executes Python/NumPy, NOT Stata or Mata.  The Stata integration suite
run_release_tests.do is the gate for the actual ado/worker/collector code.
"""
from __future__ import annotations

from dataclasses import dataclass
import itertools
import numpy as np
from numpy.testing import assert_allclose
import reselect_check as ref


def projection_selection(target, masks, eps, mu, block_size, cached):
    """Candidate blocks, preserving the first minimizer across blocks."""
    actions = np.array([ref.gls(target, mask) for mask in masks])
    delta = actions[:, ref.P:] - target[ref.P:]
    constant = np.array([target @ ref.V @ target -
                         2 * target @ ref.V @ (target - a)
                         for a in actions])
    eh = eps[:, ref.P:]
    full_projection = eh @ delta.T if cached else None
    best = np.full(len(eps), np.inf)
    chosen = np.zeros(len(eps), dtype=int)
    for left in range(0, len(masks), block_size):
        ix = np.arange(left, min(left + block_size, len(masks)))
        projection = (full_projection[:, ix] + (delta[ix] @ mu)[None, :]
                      if cached else (eh + mu) @ delta[ix].T)
        scores = projection**2 + constant[ix][None, :]
        for offset, candidate in enumerate(ix):
            improve = scores[:, offset] < best
            best[improve] = scores[improve, offset]
            chosen[improve] = candidate
    errors = eps @ target + np.sum((eh + mu) * delta[chosen], axis=1)
    return chosen, errors, best


def check_cached_blocked_scores_and_errors():
    assert np.linalg.eigvalsh(ref.V).min() > 0
    assert np.max(np.abs(ref.V[:2, 2:4])) > 0
    rng = np.random.default_rng(41026)
    draws = rng.multivariate_normal(np.zeros(len(ref.Z)), ref.V, 73)
    truth = np.r_[ref.Z[:ref.P], ref.MU]
    for target, masks in [(ref.T[0], ref.LOCAL[0]),
                          (ref.T[1], ref.LOCAL[1]), (ref.W, ref.JOINT)]:
        oracle_id, oracle_error, oracle_score = [], [], []
        for eps in draws:
            chosen, actions = ref.choose(target, masks, truth + eps)
            oracle_id.append(chosen)
            oracle_error.append(actions[chosen] @ (truth + eps) - target @ truth)
            oracle_score.append(ref.fmsc(actions[chosen], target, truth + eps))
        for cached, block in itertools.product((False, True), (1, 2, 64)):
            ids, error, score = projection_selection(target, masks, draws,
                                                      ref.MU, block, cached)
            assert np.array_equal(ids, oracle_id)
            assert_allclose(error, oracle_error, atol=1e-12, rtol=1e-12)
            assert_allclose(score, oracle_score, atol=1e-12, rtol=1e-12)


def check_candidate_ties_are_ordered():
    # Duplicate candidates deliberately give EXACT ties across separate blocks.
    masks = np.vstack([ref.LOCAL[0][0], ref.LOCAL[0][0]])
    for cached, block in itertools.product((False, True), (1, 64)):
        ids, _, _ = projection_selection(ref.T[0], masks, ref.E,
                                         ref.MU, block, cached)
        assert np.all(ids == 0)


def check_partial_cache_prefix_boundaries_with_mixed_components():
    """Global cache prefixes cross interleaved component blocks of 64."""
    rng = np.random.default_rng(730261)
    draws, candidates, dimensions = 41, 340, 4
    eh = rng.normal(size=(draws, dimensions))
    mu = rng.normal(size=dimensions)
    select_b = rng.normal(size=(candidates, dimensions))
    report_b = rng.normal(size=(candidates, dimensions))
    constants = rng.normal(size=candidates)
    component = np.arange(candidates) % 2
    baseline = rng.normal(size=draws)
    w = eh + mu
    score = (w @ select_b.T)**2 + constants
    expected = baseline.copy()
    expected_choices = []
    for group in (0, 1):
        ix = np.flatnonzero(component == group)
        choice = ix[np.argmin(score[:, ix], axis=1)]
        expected_choices.append(choice)
        expected += np.sum(w * report_b[choice], axis=1)
    for prefix in (0, 1, 63, 64, 65, 170, 339, 340):
        projection = eh @ select_b[:prefix].T
        actual = baseline.copy()
        for group in (0, 1):
            ix = np.flatnonzero(component == group)
            best = np.full(draws, np.inf)
            chosen = np.full(draws, ix[0])
            for left in range(0, len(ix), 64):
                ci = ix[left:left+64]
                cached = ci < prefix
                projected = np.empty((draws, len(ci)))
                projected[:, cached] = (projection[:, ci[cached]] +
                                         (select_b[ci[cached]] @ mu)[None, :])
                projected[:, ~cached] = w @ select_b[ci[~cached]].T
                block_score = projected**2 + constants[ci]
                for column, candidate in enumerate(ci):
                    better = block_score[:, column] < best
                    best[better] = block_score[better, column]
                    chosen[better] = candidate
            assert np.array_equal(chosen, expected_choices[group])
            actual += np.sum(w * report_b[chosen], axis=1)
        assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)


@dataclass(frozen=True)
class Result:
    start: int
    value: float
    evaluations: int
    converged: int
    exhausted: int


def budget_for_start(total, nstarts, start):
    return total // nstarts + (start <= total % nstarts)


def merge_in_original_start_order(results, objective):
    # Arrival order must not change ties.  One result per start is required.
    ordered = sorted(results, key=lambda result: result.start)
    if len({result.start for result in ordered}) != len(ordered):
        raise ValueError("duplicate start")
    direction = -1 if objective == 1 else 1
    best = ordered[0]
    for result in ordered[1:]:
        if direction * result.value > direction * best.value:
            best = result
    return (best.start, best.value, sum(r.evaluations for r in ordered),
            sum(r.converged for r in ordered), max(r.exhausted for r in ordered))


def check_budget_partition_and_worker_assignment():
    for nstarts, total in ((2, 80), (3, 83), (8, 1600), (8, 1603), (8, 8)):
        budgets = [budget_for_start(total, nstarts, s) for s in range(1, nstarts+1)]
        assert sum(budgets) == total
        assert min(budgets) >= 1
        assert max(budgets) - min(budgets) <= 1
        jobs = [(o, s) for o in range(1, 5) for s in range(1, nstarts+1)]
        for workers in (1, 2, 4, 7, 50):
            shards = [[job for index, job in enumerate(jobs)
                       if index % workers == worker] for worker in range(workers)]
            assert sorted(sum(shards, [])) == sorted(jobs)
            assert sum(len(shard) for shard in shards) == 4*nstarts


def check_out_of_order_reduction_preserves_ties_and_diagnostics():
    # Synthetic start results test the reduction contract, not optimizer quality.
    results = [Result(1, .7, 27, 0, 1), Result(2, .7, 20, 1, 0),
               Result(3, .2, 27, 0, 1)]
    for objective in range(1, 5):
        expected = merge_in_original_start_order(results, objective)
        for permutation in itertools.permutations(results):
            assert merge_in_original_start_order(permutation, objective) == expected
        assert expected[0] == (3 if objective == 1 else 1)
        assert expected[2:] == (74, 1, 1)
    try:
        merge_in_original_start_order([results[0], results[0]], 1)
    except ValueError:
        pass
    else:
        raise AssertionError("Duplicate jobs must not silently count twice")


def main():
    checks = [check_cached_blocked_scores_and_errors,
              check_candidate_ties_are_ordered,
              check_partial_cache_prefix_boundaries_with_mixed_components,
              check_budget_partition_and_worker_assignment,
              check_out_of_order_reduction_preserves_ties_and_diagnostics]
    for check in checks:
        check()
        print("PASS", check.__name__)
    print("PASS: independent algebra and scheduling checks; Stata/Mata NOT executed.")


if __name__ == "__main__":
    main()
