{smcl}
{* *! jwselect 0.4.0 23oct2026}{...}
{title:jwselect_collect -- Combine prepared focused-inference searches}

{p 4 4 2}
Requires Stata 17. Collect after the required {help jwselect_worker} jobs have
finished and released their locks. The command verifies completeness and
compatibility, combines starting points, independently checks the four winning
nuisance points, and returns marginal intervals and p-values.

{title:Syntax}

{p 8 12 2}
{cmd:jwselect_collect using} {it:"directory"} [{cmd:,} {opt progress}]

{title:Description}

{p 4 4 2}
Preparation freezes the inputs and numerical settings. Collection chooses the
minimum searched lower quantile and maxima of the other three objectives across
starts, resolving exact ties in the deterministic original order. The same
prepared independent validation draws evaluate all four quantities at all four
winning points. Exact finite-point simulation bounds use the original error
allocation. Outputs with no active nuisance adjustment use analytic normal
inference and require no numerical search jobs.

{p 4 4 2}
Independent validation currently runs serially during collection; the search
jobs provide the objective/start parallelism. The returned table retains all
point-estimate rows. Outputs not requested by
{cmd:inferoutputs()} during preparation have missing intervals and p-values.
Neither {cmd:inferoutputs()} nor job division discards the underlying joint
covariance needed by the requested target. The command does not refit raw data.

{p 4 4 2}
Collection refuses incomplete or incompatible job sets and active worker locks.
Do not combine files from different preparations. Use a new prepared directory
for changed estimates, candidates, seed, error allowances, starts, or budgets.
{cmd:progress} prints collection progress; there are no options to change frozen
inferential settings during collection.

{title:Saved results}

{p 4 4 2}
The command is r-class and leaves the current e-class estimates intact. Save
returned matrices before another r-class command overwrites them.

{synoptset 25 tabbed}{...}
{synopt:{cmd:r(table)}}estimate, se_fixed, ll, ul, pvalue for each output{p_end}
{synopt:{cmd:r(diagnostics)}}four objective diagnostics per processed output{p_end}
{synopt:{cmd:r(extrema)}}output, objective, and full nuisance means in hmap order{p_end}
{synopt:{cmd:r(output_map)}}kind and actual cohort for each output row{p_end}
{synopt:{cmd:r(jobs)}}job_id, output, objective, start{p_end}
{synopt:{cmd:r(N_jobs)}}number of numerical jobs; if zero, r(jobs) is a 1-by-4 missing placeholder{p_end}
{synopt:{cmd:r(settings)}}level, gamma, mcerr, reps, optstarts, optmaxeval, opttol, optreps, cachemb{p_end}
{synopt:{cmd:r(directory)}}prepared work directory{p_end}
{synopt:{cmd:r(inference)}}reselect{p_end}

{p 4 4 2}
Diagnostic columns: output, objective, nactive, search_reps, validation_reps,
evaluations, starts_used, converged_starts, budget_exhausted, search_value,
validated_value, validation_gap, radius_fraction. The validated value is the
independent empirical value before outward simulation correction. The final
bounds can use another objective's winning point because all four points are
checked for all four quantities.

{title:Inference limitations}

{p 4 4 2}
The finite-point independent simulation bounds do not control outer optimization
error at unvisited nuisance points. Completed jobs and mesh-converged starts do
not certify global extrema. Inspect budget/convergence diagnostics and compare
larger budgets and additional seeds. A budget-exhausted flag is not a Stata
error or a budget overrun.

{p 4 4 2}
{cmd:se_fixed} is a fixed-selection diagnostic, not a selection-adjusted normal
standard error. The returned intervals are marginal, not simultaneous bands.
The command does not return a conventional covariance matrix for the nonlinear
selected estimator. See {help jwselect} for identification and covariance-input
requirements and the numerical inference procedure.

{title:Example}

{phang2}{cmd:jwselect_collect using "C:/work/run01", progress}
{phang2}{cmd:matrix T = r(table)}
{phang2}{cmd:matrix D = r(diagnostics)}
{phang2}{cmd:matrix E = r(extrema)}
{phang2}{cmd:matrix list T}
{phang2}{cmd:matrix list D}

{title:Also see}

{p 4 4 2}
{help jwselect}, {help jwselect_worker}, {cmd:docs/jwselect_quick_user_guide.pdf}.
