{smcl}
{* *! jwselect 0.7.0 30oct2026}{...}
{title:jwselect -- Select preperiod restrictions after csdid, jwdid, or Sun-Abraham}

{p 4 4 2}
Main command for Stata 16.1 or later; the prepared-workflow commands
and release tests require Stata 17.

{title:Syntax}

{p 8 12 2}
{cmd:jwselect} [{cmd:,}
{opt select(manual|wbic|focused)} {opt wbicscope(cohort|joint)}
{opt baseline(cohort|common)} {opt anticipation(#)}
{opt focus(cohort|aggregate|both)}
{opt candidates(nested|all)}
{opt prematrix(name)} {opt candmatrix(name)} {opt jointmatrix(name)}
{opt postwindow(# #)}
{opt aggregate(equalcohort|equalcell)} {opt aggweights(name)}
{opt vce(unconditional|conditional|supplied)}
{opt restrictions(index|placebo)}
{opt level(#)} {opt inference(auto|reselect|envelope|none)}
{opt gamma(#)} {opt mcerr(#)} {opt reps(#)} {opt seed(#)}
{opt optstarts(#)} {opt optmaxeval(#)} {opt opttol(#)} {opt optreps(#)}
{opt inferoutputs(numlist)} {opt prepare(directory)} {opt cachemb(#)} {opt progress}
{opt maxcand(#)} {opt maxjoint(#)} {opt nocohort} {opt celleffects}]

{p 8 12 2}
{cmd:jwselect, bmat(B) vmat(V) cells(C) nunits(N)}
[{opt gap(#)} {opt vce(supplied)} and the selection/design options above]

{p 8 12 2}
{cmd:estat} {cmd:design}|{cmd:selection}|{cmd:scores}|{cmd:covariance}|{cmd:cells}|{cmd:inference}

{p 8 12 2}
{cmd:estat pretrend} [{cmd:,} {opt selected} {opt prematrix(P)} {opt output(#)}]

{title:Description}

{p 4 4 2}
The command obtains the unaggregated cohort-time margins and their complete
joint covariance, forms baseline treatment-effect contrasts and optional
preperiod contrasts, selects restrictions, and performs minimum-distance (MD)
adjustment. It stores complete joint fixed-selection covariance matrices.
It reads unaggregated {cmd:csdid} results or joint {cmd:jwdid estat attgt}
margins, or validated Sun-Abraham inputs from {help jwselect_sa}. Already aggregated {cmd:estat event} vectors are not supported.

{p 4 4 2}
The baseline is g-gap*(1+delta) with {cmd:baseline(cohort)} (default), or
q-gap*(1+delta) with {cmd:baseline(common)}, where q is the earliest treated
cohort in the input. For linear and matrix inputs, {cmd:anticipation(delta)}
defaults to zero. The CS-DR and nonlinear adapters retain the upstream cohort
baseline and infer delta when omitted; explicit mismatches are rejected. This differs
from jwdid's upstream {cmd:anticipation()} convention, whose default 1 omits
g-gap. Optional preperiods satisfy t<g-gap*delta and exclude the baseline.
The baseline is a maintained identifying period, never a selectable restriction.

{p 4 4 2}
Post-treatment targets use t>=g. Each displayed cohort effect is an equal-weight
mean over its included post-treatment cells. By default, all available cells
are included, so cohorts may have different exposure horizons.
{cmd:postwindow(0 2)}, for example, requires event times 0, 1, and 2 in EVERY
cohort. Missing requested support causes an error rather than a changed target.

{title:Direct linear jwdid support and covariance}

{p 4 4 2}
The initial adapter supports a balanced, regularly spaced panel, unweighted
linear fits, actual never-treated units coded gvar=0, full
{cmd:hettype(timecohort)} saturation, and clustering at {cmd:ivar()}.
The outcome and covariates must be plain numeric variables, and covariates must
be constant within individual. Explicitly generated dummy covariates are allowed.
The default reghdfe backend and {cmd:method(regress)} are supported as described
below. Supported nonlinear fits are described separately. Repeated cross-sections,
custom treatment/regression modes, weights, and other cluster designs require
separately validated matrix inputs; they are not silently accepted.

{p 4 4 2}
{cmd:vce(unconditional)} is the direct-input default and requests
{cmd:estat attgt, vce(unconditional)}. The adapter requires
{cmd:jwdid ..., method(regress) never} for this score-supported path.
{cmd:vce(conditional)} explicitly uses the usual margins covariance and also
allows jwdid's default reghdfe backend. It conditions on the observed covariate
values/distribution for the margins calculation; it does not silently provide
unconditional population-ATT inference. The command never refits the model.
Covariate averaging, when needed, occurs in jwdid/margins extraction; the
selection and MD modules then use the supplied joint estimates/covariance.

{p 4 4 2}
The adapter leaves the data unchanged. Successful jwselect estimation replaces
e() in the usual Stata manner; use {cmd:estimates store} before calling it if
you want to compare specifications. On failure, original estimation results
and the RNG state are restored. The RNG state is also restored on success.

{title:Optional Sun-Abraham adapter}

{p 4 4 2}
Fit {cmd:jwselect_sa y, ivar(id) tvar(time) gvar(first_treat)} before calling
{cmd:jwselect}. The wrapper uses optional {cmd:eventstudyinteract} and validates
the unbinned, no-covariate, never-treated design proved in draft v46_12.
It generates all observed event indicators and omits only event -1. Keep a
stored source fit and restore it before each selection specification.
See {help jwselect_sa} for installation, supported data, and saved results.

{p 4 4 2}
The full unit-clustered regression covariance is used with {cmd:vce(supplied)}.
Linear rebasing via {cmd:baseline(common)} or {cmd:anticipation()} is supported.
All reporting weights remain fixed: population IW aggregation with estimated
cohort-share uncertainty is not implemented. Importing raw, unverified
{cmd:eventstudyinteract} or {cmd:reghdfe} results is not supported.

{title:Optional CS-DR adapter}

{p 4 4 2}
Install {cmd:csdid} separately using {browse "https://psantanna.com/csdid/":the authors' installation guide}.
The adapter targets the 2.0.1 API and requires version 2.0.1 or later in the 2.x
line; older SSC 1.x releases are rejected. The pinned 2.0.1 source commit is
{cmd:f98bcea0691b1670c409bb7e3cc81d60705b4f84}; the README gives its installation
command. Record {cmd:csdid version} for reproducibility. No dependency is installed by
{cmd:jwselect}; the matrix interface needs no upstream estimator.

{p 4 4 2}
Run {cmd:csdid} with {cmd:method(dr) nevertreated base_period(universal) analytical pointwise}.
The supported sample is an unweighted balanced panel with positive integer,
unit-spaced calendar time, actual never-treated controls coded zero, and plain
numeric covariates constant within unit. Omit extra clustering or use
{cmd:cluster(ivar)}; higher-level clusters are unsupported. Keep the original
estimation data. Call {cmd:jwselect} on unaggregated ATT(g,t) results, before
posting an aggregation. Bootstrap covariance, varying bases, not-yet-treated
controls, and unsupported samples are rejected.

{p 4 4 2}
{cmd:vce(supplied)} is the only CS-DR covariance convention and is the default
for that source. The adapter retains the full analytical influence-function
covariance, including cross-cohort and cross-period terms. It reconstructs
structural-zero reference cells from the upstream cell metadata.
{cmd:baseline(cohort)} is required; an omitted {cmd:anticipation()} is inferred
from {cmd:csdid}'s value, and an explicit value must match it.

{title:Optional nonlinear jwdid adapter}

{p 4 4 2}
{cmd:method(logit)} and {cmd:method(poisson)} are supported under the balanced,
unweighted, never-treated, unit-clustered and time-invariant-covariate contract
above. Other nonlinear methods, including {cmd:ppmlhdfe}, are not supported.
{cmd:restrictions(index)} selects the full lead coefficient block: the intercept
and every covariate interaction at a cohort/date. The ATT remains on the outcome
scale. Unit scores provide the full joint covariance between ATT and index
coefficients, retaining shared nuisance and cross-cohort/date terms.
{cmd:vce(unconditional)}, the default, adds empirical target-covariate-distribution
uncertainty and its cross terms; {cmd:vce(conditional)} conditions on the observed
covariates. Centered interactions and upstream {cmd:xasis} are supported; index
coefficients are reported in the original covariate basis.

{p 4 4 2}
Each complete date block is one selection unit. With q covariates it has q+1
identified coefficients. Wald-BIC penalizes the restriction rank, not the number
of selected dates. Generated candidate families enumerate dates; {cmd:prematrix()}
selects complete blocks. Custom candidate masks have {cmd:e(K)} coefficient
columns in {cmd:e(hmap)} order and must assign the same bit to every component
of a date. Split, unidentified, or singular blocks are rejected.

{p 4 4 2}
For backward compatibility, omitted {cmd:restrictions()} means
{cmd:restrictions(placebo)}: joint {cmd:estat attgt} outcome-scale ATT and average
placebo margins supply the existing inputs. The new option is only accepted
after supported nonlinear jwdid logit or Poisson fits.

{p 4 4 2}
Retain {cmd:baseline(cohort)} and the upstream normalization. If
{cmd:anticipation()} is omitted, the adapter uses the recorded {cmd:jwdid}
anticipation value minus one. An explicit incompatible value, or
{cmd:baseline(common)}, is rejected. The original linear/matrix rebasing
behavior is unchanged.

{title:Moment interpretation for CS-DR and nonlinear inputs}

{p 4 4 2}
The maintained estimate d is an outcome-scale ATT; optional h values are zero
average outcome-scale placebo effects, or full conditional-index lead blocks
under {cmd:restrictions(index)}. For a fixed candidate, the adjustment
is d - Cov(d,h) inv(Var(h)) h, using the same joint covariance throughout.
This is a Gaussian covariance-weighted moment adjustment. Its Wald--J identity
belongs to the canonical quadratic MD problem with unrestricted post effects
and h=0. It does not assert equality with a constrained nonlinear refit or
CS-DR re-estimation at alternative base periods.

{p 4 4 2}
With covariate interactions, zero average nonlinear placebo effects can be
weaker than setting every index-pretrend coefficient to zero. The index route
implements the v48_3 full-block target-level MD construction; the separate
coefficient-MD projection followed by nonlinear prediction is not implemented. Identification,
joint asymptotic normality and consistent covariance are required. CS-DR
robustness conditions apply to the baseline and each placebo moment used;
the backend name alone does not establish those conditions. The MD adjustment
can leave natural outcome bounds; results are not clipped. Fixed
{cmd:jwselect} aggregation weights need not reproduce the upstream commands'
population-share aggregates. Upstream confidence intervals do not adjust for
{cmd:jwselect}'s selection step.

{title:Selection}

{p 4 4 2}
{cmd:select(manual)} uses exactly the optional dates in {cmd:prematrix(P)}.
P has two columns: cohort, calendar time. Do not list baseline dates.
Omitting P means baseline-only estimation for every cohort.

{p 4 4 2}
{cmd:select(wbic)} is the default. With {cmd:wbicscope(cohort)} (default),
select a candidate within each cohort by minimizing Wald - r*ln(N), where r
is the number of optional restrictions and N is the independent-unit count.
The final cohort/cell/aggregate estimates use JOINT MD with the union of the
selected restrictions. Thus another cohort's restriction can affect a cohort
estimate through the joint covariance.

{p 4 4 2}
{cmd:select(wbic) wbicscope(joint)} minimizes the same criterion over complete
joint configurations, using the full covariance across cohorts and dates.
The joint Wald statistic equals the covariance-compatible MD J statistic.
The selected configuration is then used in joint MD for every output.
Joint WBIC can select different restrictions than cohort-wise WBIC.
{cmd:wbicscope()} is a WBIC selection option; it does not change focused FMSC.

{p 4 4 2}
{cmd:select(focused)} minimizes FMSC = v0 - 2*gain + adjustment^2 in the
sampling-covariance scale. This is the draft's asymptotic-covariance score
divided by N; it gives the same minimizer. The squared bias correction is
untruncated. Exact focused-score ties favor fewer restrictions, then the first stored candidate among equal sizes. Generated and custom
local families put baseline-only first; a supplied jointmatrix() retains its
first-occurrence row order and prepends baseline-only only if absent.

{p 4 4 2}
{cmd:focus(cohort)} selects for each cohort's displayed mean using that cohort's
candidate restrictions, applies that selection to its cells, and reports the
weighted aggregate of those selected cell estimates. {cmd:focus(aggregate)}
selects a JOINT restriction configuration for the aggregate and reports cohort
estimates under that common selected joint model. These cohort rows are not
cohort-focused estimates. {cmd:focus(both)} (default) reports cohort-focused
rows, their aggregate ({cmd:agg_cohort}), and a separately optimized aggregate
({cmd:agg_focused}). The two aggregate estimators have the same fixed target.

{p 4 4 2}
{cmd:candidates(nested)} starts with the empty set and adds eligible dates in
reverse chronological order, one at a time within each cohort.
{cmd:candidates(all)} enumerates all subsets. Neither family is chosen by
preliminary significance tests. All candidate families and weights should be
fixed before inspecting outcomes.

{p 4 4 2}
{cmd:candmatrix(CM)} supplies custom local families. CM's first column is the
actual cohort and its remaining K columns are binary masks in the order of
{cmd:e(hmap)}. A local candidate can select only its own cohort's columns.
Supply at least one row for every cohort. The empty candidate is inserted first;
duplicate nonempty candidates are rejected. Run baseline-only manual estimation
and inspect {cmd:estat design} to obtain the map, then restore the upstream fit or reuse
matrix inputs for the final call.

{p 4 4 2}
Joint WBIC and aggregate-focused search use the Cartesian product of local families unless
{cmd:jointmatrix(JM)} supplies a prespecified binary matrix with K columns in
hmap order. Each row of JM is a complete joint configuration; it need not be
a Cartesian-product member. The empty row is prepended only if absent; duplicates
are coalesced in first-occurrence order. Focused-score ties favor fewer restrictions, then stored order; Wald-BIC ties use stored order. {cmd:maxcand(4096)} limits each local family and
{cmd:maxjoint(4096)} limits the joint family. Exceeding a limit causes an error;
there is no truncation or greedy approximation. For K_g eligible optional dates per cohort, the joint
Cartesian family has product(K_g+1) nested configurations or 2^sum(K_g)
all-subset configurations. The joint cap can bind even when every local
family is below maxcand(). Explicit candidate loading storage can be substantial.

{title:Joint pretrend tests}

{p 4 4 2}
{cmd:estat pretrend} jointly tests all eligible optional preperiod contrasts
(or outcome-scale placebo moments / full index blocks for the adapters) using their full covariance. {cmd:estat pretrend, prematrix(P)} tests
a prespecified subset; P has two columns, cohort and calendar time, and must
contain only eligible optional dates. The maintained baseline is never tested.

{p 4 4 2}
{cmd:estat pretrend, selected} tests the restriction mask for a reported output.
For focused results, {cmd:output(#)} is REQUIRED with {cmd:selected}; inspect
{cmd:e(output_map)} to identify the row. For manual/WBIC results it defaults
to output 1 because all outputs share the same joint restriction set.
{cmd:output()} is allowed only with {cmd:selected}. The options {cmd:selected}
and {cmd:prematrix()} are mutually exclusive. For {cmd:agg_cohort}, the
selected mask is the union of the cohort-selected restrictions. The reported
J is the compatible joint-MD criterion for that tested set, not a residual
criterion evaluated at the separately adjusted focused estimate.

{p 4 4 2}
For restriction contrasts h and their sampling covariance V_h, the command
computes Wald = h' inv(V_h) h, and reports J equal to Wald using the
covariance-weighted MD identity. J is not independently computed by refitting
or minimizing a second model. The equality requires the compatible MD model, with unrestricted postperiod
effects and the same covariance construction. The test checks that the
chosen restriction covariance is positive definite; returning J through the
identity does not independently certify the identification or full-covariance
assumptions. With r prespecified restrictions,
the reference distribution is chi-square(r) under the maintained assumptions.
A selected-set p-value is a diagnostic, not automatically valid after
outcome-dependent selection; the same caution applies to a data-chosen custom
P. This command neither tests the baseline identifying assumption nor supplies
postselection inference for treatment effects.

{p 4 4 2}
Saved scalars are {cmd:r(chi2)}, {cmd:r(J)}, {cmd:r(df)}, {cmd:r(p)},
{cmd:r(postselection)} (1 for {cmd:selected}, otherwise 0), and {cmd:r(output)}
(selected output row, or 0 for all/custom).
{cmd:r(test_set)} is all, selected, or custom. {cmd:r(table)} has columns
{cmd:df Wald J pvalue}. {cmd:r(mask)} records the tested mask in hmap order
when K>0; {cmd:r(preperiods)} lists tested cohort/date pairs when nonempty.
For an empty test, df=0 and Wald=J=0, with p missing. The command preserves
active estimation results, data, and RNG state.

{title:Aggregation weights}

{p 4 4 2}
{cmd:aggregate(equalcohort)} is the default: every cohort gets weight 1/G and
its included post cells share that weight equally. {cmd:aggregate(equalcell)}
weights all included post cells equally. Both use fixed weights conditional
on the specified target design.

{p 4 4 2}
{cmd:aggweights(W)} overrides these rules. W has columns cohort, time, weight;
it must list every included target cell exactly once, including explicit zeros.
Weights must be nonnegative and sum to one. They are never renormalized.
With custom within-cohort weights, the aggregate is a weighted sum of cell
effects and may not be a weighted sum of the displayed uniform cohort means.
Estimated population-share weights require extra influence-function terms and
are not implemented; the command treats supplied weights as fixed. Signed
contrasts are outside this initial averaging interface.

{title:Inference}

{p 4 4 2}
Manual selection reports normal inference under the fixed maintained
restrictions. Wald-BIC reports pointwise oracle normal inference; this is not
uniform inference under local misspecification.

{p 4 4 2}
For focused selection, {cmd:inference(auto)} now resolves to
{cmd:inference(reselect)}. At every trial nuisance value, common joint Gaussian
draws rerun the actual focused candidate rule. The selected estimation error
includes the nuisance mean. Selection and reporting targets are distinguished:
aggregate-focused cohort rows use the aggregate's selection rule, and custom
weighted aggregates of cohort-focused estimates retain each cohort's original
selection rule. Aggregates are formed within each draw before taking quantiles,
preserving cross-cohort covariance.

{p 4 4 2}
The two-step procedure builds an ellipsoid for the precontrast nuisance mean
and uses native Mata multistart projected pattern search over its standardized
unit ball. This is not NOMAD and needs no external optimizer. Separate searches
minimize the lower error quantile, maximize the upper quantile, and maximize
the two zero-null tail probabilities. {cmd:optstarts(4)} sets the number of
starts per objective; {cmd:optmaxeval(400)} sets its total evaluation budget;
{cmd:opttol(.02)} sets the stopping mesh in unit-ball coordinates;
{cmd:optreps(1999)} sets the common search-draw count.
The defaults are starting numerical budgets, not accuracy guarantees.

{p 4 4 2}
{cmd:gamma(.01)} allocates error to the nuisance ellipsoid.
{cmd:mcerr(.001)} allocates error to an independent simulation check of the
four discovered points. {cmd:reps(19999)} sets the independent checking-draw
count and {cmd:seed(12345)} controls the RNG. Exact binomial order-statistic
quantile bounds and one-sided Clopper-Pearson tail bounds use Bonferroni
allocation over 16 finite-point bounds per output. The remaining nominal
error is 1-level/100-gamma-mcerr. Infeasible simulation sizes cause an error;
increase {cmd:reps()}.

{p 4 4 2}
A numerical search cannot certify the global nuisance extrema. Independent
simulation controls error at discovered points only; it does not cover
unvisited nuisance points. Thus this is a numerical implementation of the
paper's theoretical construction, without a certified worst-case coverage
guarantee. {cmd:estat inference} reports evaluation counts, starts, stopping
status, searched values, independent bounds, and discovered-point locations.
Budget exhaustion returns the best point found with a warning. Increasing
starts, evaluations and draws and reducing the mesh tolerance checks numerical
stability. Even successful mesh stopping is not proof of global optimization.

{p 4 4 2}
{cmd:inference(envelope)} explicitly retains the earlier all-candidate analytic
bound with its DKW simulation correction. It ignores whether a candidate wins
selection within a draw and can be substantially more conservative.
{cmd:inference(reselect)} and {cmd:inference(envelope)} require focused selection.
{cmd:inference(none)} skips CIs and p-values. Neither approach guarantees
shorter intervals than baseline. Intervals are marginal, not simultaneous
bands, and rely on the joint Gaussian/local asymptotic approximation and the
supplied covariance; they are not finite-sample guarantees.

{p 4 4 2}
All displayed SEs are labeled {cmd:se_fixed}. Under focused selection they are
fixed-action diagnostics, not standard errors for a normal approximation to
the selected estimator. Accordingly, focused results post e(b) without e(V).
Use e(ci) and e(pvalue). e(V_fixed) and e(V_cells) preserve the joint
fixed-selection covariance for diagnostics; do not pass them to conventional
{cmd:lincom}, {cmd:test}, or {cmd:margins} as selection-adjusted covariance.
No single covariance matrix characterizes the focused estimator's nonnormal
selection distribution. Observation-level {cmd:predict} is deliberately undefined.

{title:Requested outputs, caching, and prepared jobs}

{p 4 4 2}
{cmd:inferoutputs(numlist)} identifies output ROW NUMBERS in e(output_map),
not cohort dates or event times. It requires focused selection with
{cmd:inference(none)}, {cmd:auto}, or {cmd:reselect}. Without it, inference is
requested for every output. Unrequested outputs retain point estimates and
fixed-selection covariance but have missing intervals and p-values. Use this
option to avoid an unused aggregate's numerical inference. {cmd:nocohort}
changes display only. {cmd:focus(cohort)} also avoids asking for separate
aggregate-focused point selection when that estimator is not needed.

{p 4 4 2}
{cmd:cachemb(64)} sets the search-projection cache allowance in decimal MB per
active calculation. It caches as many initial candidate-projection columns
as fit and streams the remainder in original order, avoiding an all-or-nothing
threshold. {cmd:cachemb(0)} streams all candidates. This is not a total RAM limit. Independent validation is blocked rather than keeping a complete large
validation projection cache. Account for all concurrent workers when choosing
memory budgets. {cmd:progress} requests numerical progress messages.

{p 4 4 2}
{cmd:prepare("newdir")} requires {cmd:select(focused) inference(none)} and a
NEW directory whose parent exists. It saves point-selection inputs, full joint
covariance, candidate definitions, settings, common search draws, independent
validation draws, and fixed starting points. It does not execute numerical
search. Each active requested output creates four objectives times the number
of starts. Outputs needing no active nuisance search are handled analytically.

{p 4 4 2}
Run {help jwselect_worker} in one or more Stata 17 sessions, then
{help jwselect_collect}. For example:

{phang2}{cmd:jwselect, select(focused) focus(cohort) postwindow(0 0) inference(none) inferoutputs(1) prepare("C:/work/newrun") optstarts(8) optmaxeval(1600) cachemb(64)}
{phang2}{cmd:jwselect_worker using "C:/work/newrun", worker(1) workers(1) cachemb(64) progress}
{phang2}{cmd:jwselect_collect using "C:/work/newrun", progress}

{p 4 4 2}
For four concurrent sessions, use worker indices 1 through 4 with workers(4).
Each may be further filtered by output, objective, starting point, or original
job IDs through {cmd:jobids()}; filters are AND-combined with the original-ID
partition. Batching assigned IDs for a target reuses its loaded cache/context. The total
search budget remains PER OBJECTIVE across starts: 1,600 evaluations and eight
starts give 200 evaluations per start. More sessions do not increase it.
Prepared configurations and per-job locks reject incompatible or overlapping
work. Completed compatible jobs are reusable. Use a new prepared directory for
changed settings, seeds, inputs, or release versions.

{p 4 4 2}
The collector returns r(table), r(diagnostics), r(extrema), r(output_map),
r(jobs), r(N_jobs), and r(settings). It does not replace e(V) with a selected-
estimator covariance. Serial/divided equivalence requires identical inputs,
draws, starts, budgets, tie order, and cache setting. Caching versus streaming
uses the same formulas, but rounding near ties may change a search trajectory;
check differences using numerical tolerances and diagnostics.

{title:Matrix interface}

{p 4 4 2}
B is a row or column vector of unaggregated cohort-time estimates; V is its
sampling covariance, NOT N times its sampling covariance. C has corresponding
rows (cohort,time). Keep normalization rows, even if their estimates and
covariance rows are zero. Raw V may be positive semidefinite; each candidate's
restriction covariance must be nonsingular. Invalid or unidentified
restrictions are rejected, not handled by a silent generalized inverse.
Rows may arrive in any order. Only positive treated-cohort codes belong in C.
{cmd:gap(1)} defaults to one; times must align to this period scale.
{cmd:nunits(N)} declares the independent-unit count. The user is responsible
for the estimand and covariance convention of externally supplied matrices.

{title:Saved results}

{p 4 4 2}
Key matrices: e(b), e(V) for manual/WBIC only, e(V_fixed), e(table), e(ci),
e(pvalue), e(output_map), e(baselines), e(postmap), e(hmap), e(selection),
e(cell_selection), e(b_cells), e(V_cells), e(weights), e(scores),
e(candidates), e(joint_candidates), e(z), e(V_z), e(loading), e(targets).
Numerical reselection also stores e(inference_diagnostics) and
e(inference_extrema).
If K=0, hmap/selection/cell_selection/joint_candidates are absent.
e(joint_candidates) is present for joint WBIC and aggregate-focused searches
with K>0.

{p 4 4 2}
Canonical z orders post targets first, optional restriction estimates second.
For index mode these are outcome-scale ATTs followed by raw index coefficients.
e(loading)*e(z)' reproduces e(b)'; e(loading)*e(V_z)*e(loading)' reproduces
e(V_fixed). e(output_map) columns are kind and actual cohort. Kinds: 1 cohort
(manual/WBIC/cohort focus), 2 aggregate of cohort-focused cells, 3 aggregate
focus, 4 manual/WBIC aggregate, 5 cohort under aggregate-focused joint fit.
Cohort coefficient names are cohort1, cohort2, ...; use output_map for dates.

{p 4 4 2}
e(scores) columns: scope (1 local, 2 joint), cohort (0 joint), candidate,
nrestr, Wald, WBIC, estimate, se_fixed, adjustment, gain, FMSC.
Joint WBIC rows have scope=2 and cohort=0; selection minimizes the WBIC
column, while focused selection minimizes FMSC. The other columns are diagnostics.
Local score estimates use local restrictions even when the final WBIC/manual
fit uses the joint union. e(candidates) columns: cohort, local candidate ID,
then K mask bits. e(selection) records all restrictions used by each output.

{p 4 4 2}
Scalars include N_units, N_observations, N_covariates, N_clusters, G, K, P,
gap, anticipation, upstream_anticipation, level, gamma, mcerr, reps, seed,
maxcand, maxjoint, optstarts, optmaxeval, opttol, optreps, cachemb. Matrix input leaves
unavailable data counts missing.
Macros record baseline, selection_method, wbic_scope, focus, candidate_type, aggregate,
covariance, inference, postwindow, sourcecmd, cmdline, inferoutputs,
and prepared_directory. Source metadata are e(source_estimator),
e(source_version), e(source_method), e(source_model), and e(restriction_scale).
The source estimator is matrix, csdid, jwdid, or eventstudyinteract. e(source_version) is empty
when the source does not expose a verified version, including jwdid; retain
{cmd:which jwdid} output separately. e(version) and e(build) identify the
installed release.

{p 4 4 2}
For each output, e(inference_diagnostics) has four rows, one per
objective (1 lower quantile; 2 upper quantile; 3 left tail; 4 right tail).
Columns are output, objective, nactive, search_reps, validation_reps,
evaluations, starts_used, converged_starts, budget_exhausted, search_value,
validated_value, validation_gap, radius_fraction. validated_value is the
independent empirical objective before Monte Carlo bound correction;
validation_gap is its absolute difference from search_value. The radius
fraction is the unit-ball norm. e(inference_extrema) contains output, objective, then K nuisance
means in hmap order; inactive means equal the observed restriction estimates.
Final bounds use extrema over all four checked points, so an output endpoint
need not use the point from the matching objective row. Normal-inference
outputs with no active adjustment have zero counts and missing objectives.

{title:Conditional-index stored results}

{p 4 4 2}
{cmd:e(K)} counts scalar index restrictions; {cmd:e(K_blocks)} counts date blocks.
{cmd:e(blockmap)} contains cohort, time, and rank. {cmd:e(index_b)},
{cmd:e(V_index)}, and {cmd:e(Cov_att_index)} contain lead coefficients, their
covariance, and raw-ATT/index covariance. {cmd:e(index_hmap)} maps components to
cohort/date; {cmd:e(upstream_cells)} maps raw ATT rows. {cmd:e(index_terms)} records
source coefficient positions and components (0 intercept, 1..q covariates in
{cmd:e(index_covariates)} order). For centered fits the intercept position is an
anchor: the reported raw intercept also subtracts cohort means times slopes.
{cmd:e(index_basis)} is raw_covariates. {cmd:estat design} displays these maps.
{cmd:e(index_covariance_factor)} is the validated common finite-sample factor;
{cmd:e(index_covariance_error)} and {cmd:e(index_margins_error)} record relative
reconstruction errors. {cmd:e(upstream_vce)} is conditional for the source table;
adjusted inference uses the requested {cmd:e(covariance)} convention.
With no eligible leads, index matrices are absent and both counts are zero.

{title:Examples}

{phang2}{cmd:csdid y x1 x2, ivar(id) time(year) gvar(first_treat) method(dr) nevertreated base_period(universal) analytical pointwise cluster(id)}
{phang2}{cmd:estimates store CS_SOURCE}
{phang2}{cmd:jwselect, select(wbic) postwindow(0 2)}
{phang2}{cmd:estimates restore CS_SOURCE}
{phang2}{cmd:jwselect, select(manual) postwindow(0 2)}

{phang2}{cmd:jwdid y x1 x2, ivar(id) tvar(year) gvar(first_treat) never method(logit) cluster(id)}
{phang2}{cmd:estimates store NL_SOURCE}
{phang2}{cmd:jwselect, select(wbic)}
{phang2}{cmd:estimates restore NL_SOURCE}
{phang2}{cmd:jwselect, select(manual)}

{p 4 4 2}
For an exponential conditional mean use {cmd:method(poisson)} in the source fit.
Self-contained examples are {cmd:examples/csdid_demo.do} and
{cmd:examples/nonlinear_demo.do}. Restore the stored source fit before each
new direct specification. The following example retains the linear workflow.


{phang2}{cmd:jwdid y x1 x2, ivar(id) tvar(year) gvar(first_treat) never method(regress)}
{phang2}{cmd:estimates store original}
{phang2}{cmd:jwselect, select(wbic) wbicscope(joint) candidates(nested)}
{phang2}{cmd:estat pretrend}
{phang2}{cmd:estat pretrend, selected}
{phang2}{cmd:estimates restore original}
{phang2}{cmd:jwselect, select(focused) focus(both) postwindow(0 2)}
{phang2}{cmd:estat selection}
{phang2}{cmd:estat pretrend, selected output(1)}
{phang2}{cmd:estat covariance}
{phang2}{cmd:estat inference}
{phang2}{cmd:jwselect, nocohort}

{p 4 4 2}
For a numerical stability check, restore the same upstream estimates and use
{cmd:select(focused) optstarts(8) optmaxeval(1600) opttol(.01) optreps(9999) reps(49999)}.
A separate {cmd:inference(envelope)} run compares the previous conservative
construction while leaving point estimates unchanged. Large all-subsets
families can make these checks expensive.

{p 4 4 2}
For reproducible matrix-input examples, see {cmd:examples/matrix_demo.do}
and {cmd:examples/joint_demo.do}. Both examples are self-contained.

{p 4 4 2}
Developer tests are supplied separately in the full source archive, not this
clean distribution. From that source archive, the joint public-command test is
{cmd:do tests/joint040_smoke.do}, also included in the release driver.
Run {cmd:do run_release_tests.do} from the package directory in a fresh Stata 17
session. The test clears data and estimates; save your work first.
Require the final release PASS. No test installs or updates packages automatically.
