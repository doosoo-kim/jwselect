# Sun–Abraham Monte Carlo for jwselect 0.7.0

`jwselect_sa_monte_carlo.do` measures bias, RMSE, and confidence interval coverage for the Sun–Abraham adapter, including after preperiod selection. It refits `jwselect_sa` in every replication and uses the full joint covariance returned by the wrapper. The theory and interpretation follow draft **v46_14**. This is a new simulation design for the adapter, not a reproduction of the paper's simulation tables.

Run it in a fresh Stata 16.1+ session: the do-file clears the current data and estimation results. Install the package and its documented optional dependencies (`eventstudyinteract`, `reghdfe`, `ftools`, and `avar`) beforehand. The file installs nothing and does not change your installed package. The supported `eventstudyinteract` source is checked by `jwselect_sa`.

From the unpacked package directory:

```stata
do examples/jwselect_sa_monte_carlo.do smoke 2 600 "." "."
do examples/jwselect_sa_monte_carlo.do pilot 100 600 "." "."
do examples/jwselect_sa_monte_carlo.do full 1000 1000 "." "."
```

Arguments are `mode`, replications **per DGP**, number of units, package directory, and existing output parent. Paths can contain spaces if quoted. Use an absolute path to the do-file and supply an absolute package directory when running from elsewhere. N must be at least 100 and divisible by 10. The default with no arguments is `full`, 1,000 replications per DGP, N=1,000; this is a substantial computation. Start with `smoke` to check the installation.

| Profile | Outer reps per DGP | N | Starts per objective | Total evaluations per objective | Stopping mesh | Search draws | Independent checking draws |
|---|---:|---:|---:|---:|---:|---:|---:|
| smoke | 2 | 600 | 2 | 200 | .05 | 199 | 999 |
| pilot | 100 | 600 | 4 | 800 | .02 | 999 | 9,999 |
| full | 1,000 | 1,000 | 8 | 6,400 | .01 | 9,999 | 49,999 |

Command-line values override the outer replication count and N. Edit the numerical settings in the profile block to compare different search and simulation budgets. The candidate family, target, and DGP formulas are the same across profiles. Reuse N and outer replication numbers to compare budgets on identical simulated panels.

## Data-generating process and target

There are eight periods and fixed cohort counts: 40% never treated, 40% first treated at period 4, and 20% first treated at period 6. Unit effects vary across cohorts. Errors follow a stationary Gaussian AR(1) process with correlation .5 and marginal standard deviation .8. The treatment effect in cohort g at event time e is

`tau(g,e) = .5 + .15 * 1(g=6) + .08 * e`.

The target is the equally weighted average of the two cohorts' mean treatment effects over event times 0, 1, and 2:

`(.58 + .73) / 2 = .655`.

Weights are prespecified at one half per cohort, and one third per event time within cohort. This is a conditional-design experiment with fixed cohort counts and fixed aggregation weights. It does not assess population interaction-weighted inference with estimated cohort shares.

The three DGPs add an untreated-outcome shift to treated cohorts in periods 1 and 2 only:

| Design | Early-period shift | Interpretation |
|---|---:|---|
| valid | 0 | All optional preperiod restrictions hold. |
| local | `3/sqrt(N)` | A fixed local misspecification vector as N increases. |
| fixed | .5 | A fixed violation; useful as a stress test outside the local-sequence coverage theorem. |

The omitted baselines, period 3 for cohort 4 and period 5 for cohort 6, remain valid in every DGP. The postperiod target remains .655. Optional contrasts are `(4,1), (4,2), (6,1), (6,2), (6,3), (6,4)`. All 64 subsets are available for joint WBIC and aggregate-focused selection. The no-restriction baseline is included.

All five methods use the same simulated panel within a replication. The same underlying innovations also pair the three DGPs. The data seed is `614729 + replication`; inference has a separate seed, `371829 + replication`. The generator is explicitly `mt64`. Changing inference budgets does not change the simulated data.

## Methods and interval interpretation

| Output method | Point estimate | Interval |
|---|---|---|
| baseline | Prespecified no-optional-restriction estimate | Ordinary normal interval returned by `jwselect`. |
| all_pre | Prespecified all-six-restriction MD estimate | Ordinary normal interval; intentionally misspecified in the invalid designs. |
| wbic_joint_normal | Joint Wald–BIC selected estimate | Pointwise-oracle normal interval. Its local-misspecification coverage is an empirical benchmark, not a uniform coverage claim. |
| focused_normal | Aggregate-focused selected estimate | Deliberately naive normal interval from the fixed-action SE. |
| focused_reselect | Exactly the same estimate as `focused_normal` | `inference(reselect)` interval extracted from `e(ci)`, with the rule rerun in every Gaussian draw. |

Only the equal-cohort aggregate is evaluated for coverage; `inferoutputs(3)` avoids computing unused cohort intervals. The script verifies the output map before extracting this row. `jwselect 0.7.0` supports reselection inference for focused selection only; it has no WBIC reselection interface.

Nominal marginal coverage is 95%, with `gamma(.01)` for the nuisance region and `mcerr(.001)` for checking the four discovered points. Independent checking draws provide finite-point Monte Carlo bounds. They do not certify the nuisance extrema over the entire ellipsoid. Increasing numerical budgets and checking stability is necessary when evaluating coverage; even successful mesh stopping is not proof of global optimality. Report the numerical settings alongside coverage results.

## Saved output and diagnostics

Every run creates a unique `sa_mc_<date>_<time>_<suffix>` directory:

- `settings.txt` and `run.log`: configuration, command versions, progress, and totals.
- `sa_mc_draws.dta` / `.csv`: every attempted DGP–replication–method row, including failures, true target, point estimate, fixed-action SE, CI endpoints, selection mask, seeds, return codes, and optimizer summaries.
- `inference_diagnostics.dta` / `.csv`: four objective rows for every successful focused reselection calculation, including active dimension, counts, convergence, exhaustion, search/check values, validation gaps, and nuisance-region locations as radius fractions.
- `sa_mc_summary.dta` / `.csv`: bias, RMSE, coverage, their Monte Carlo standard errors, exact binomial Monte Carlo coverage limits, mean interval width, mean number of selected restrictions, failure denominators, and optimizer-flag counts.

`source_rc`, `point_rc`, and `inference_rc` are separate. A successful focused point estimate is retained even if adjusted inference fails. No failed replication is retried or silently discarded. `selected_mask` uses bits in the six-contrast order above; zero is the baseline, and 63 means all six restrictions.

Bias and RMSE use `n_estimate` successful point estimates. Coverage and mean width use `n_interval` valid returned intervals; the failure count and total attempts are shown separately. `covered_per_attempt` counts unavailable intervals as failure to deliver coverage and divides by all attempted replications. This is an availability-adjusted rate, not ordinary conditional-on-returned-interval coverage. If any failures occur, report both rates and investigate their return codes before interpreting performance.

The standard errors are `sd(error)/sqrt(n_estimate)` for bias, the delta-method `sd(error^2)/(2*RMSE*sqrt(n_estimate))` for RMSE, and `sqrt(p*(1-p)/n_interval)` for coverage. Exact binomial 95% Monte Carlo limits accompany coverage because the plug-in standard error is zero when a tiny simulation happens to produce coverage 0 or 1. Neither smoke nor a small pilot is evidence of nominal coverage. At true coverage .95, 1,000 replications have an approximate standard error of .0069, or .69 percentage points.

`n_budget_flag` counts intervals with at least one exhausted objective. `n_mesh_flag` counts intervals with at least one unconverged attempted search start. Such intervals are retained in the main summary, so these flags cannot improve reported coverage through selective exclusion. Increase evaluations/starts and repeat on the same panels if flags are frequent. The detailed objective-specific validation gaps are more useful than the scalar maximum: quantile gaps and tail-probability gaps have different units. Search convergence and small gaps still do not provide a global optimization certificate.

This simulation tests implementation and finite-sample performance under the stated DGPs. It does not prove the paper's asymptotic guarantee, general validity of arbitrary data-dependent preperiod selection, or coverage outside the supported Sun–Abraham adapter design.
