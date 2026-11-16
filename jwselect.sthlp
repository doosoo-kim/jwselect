{smcl}
{* *! jwselect 0.9.1 16nov2026}{...}
{title:jwselect -- Select preperiod restrictions after csdid, jwdid, or Sun-Abraham}

{p 4 4 2}
The main command requires Stata 16.1 or later. Prepared workers and collection require Stata 17 or later.

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
{opt level(#)} {opt inference(auto|wbic|reselect|envelope|none)}
{opt gamma(#)} {opt mcerr(#)} {opt reps(#)} {opt seed(#)}
{opt wbicmethod(bootstrap|mc)} {opt wbicengine(integrated|direct)}
{opt mcfallback(terminal|error)} {opt mcbeta(#)} {opt mcsteps(#)} {opt mcrefine}
{opt mccells(#)} {opt mcvalidate(#)} {opt mcradius(#)} {opt mcpmax(#)}
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

{p 8 12 2}
{cmd:estat compare} [{cmd:,} {opt cells} {opt output(#)}]

{title:Description}

{p 4 4 2}
The command obtains the unaggregated cohort-time margins and their complete
joint covariance, forms baseline treatment-effect contrasts and optional
preperiod contrasts, selects restrictions, and performs minimum-distance (MD)
adjustment. It stores complete joint covariance matrices calculated by treating
the selected restrictions as fixed.
It reads unaggregated {cmd:csdid} results or joint {cmd:jwdid estat attgt}
margins, or validated Sun-Abraham inputs from {help jwselect_sa}.
The scalar sources {help jwselect_rcs} and {help jwselect_sa_iw} instead supply
one prespecified population target, optional preperiod moments, and their
complete joint covariance. Already aggregated {cmd:estat event} vectors are not supported.

{p 4 4 2}
This release documents the manuscript dated November 16, 2026. Retain the
maintained baseline estimate as a reference and report MD estimates alongside
it. Asymptotic efficiency improvement does not guarantee finite-sample
improvement for an unknown data-generating process. {cmd:estat compare}
provides a matched baseline, own-cohort MD, and joint MD comparison.

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
These design options apply to ordinary cohort-time sources. After a scalar
source, set the target, baseline, population and event window in that source
command; {cmd:jwselect} cannot change them. Scalar sources report only one
population estimate, and focused selection uses {cmd:focus(aggregate)}.

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
an unbinned design with never-treated controls and no covariates.
It generates all observed event indicators and omits only event -1. Keep a
stored source fit and restore it before each selection specification.
See {help jwselect_sa} for installation, supported data, and saved results.

{p 4 4 2}
The full unit-clustered regression covariance is used with {cmd:vce(supplied)}.
Linear rebasing via {cmd:baseline(common)} or {cmd:anticipation()} is supported.
This direct route uses fixed reporting weights. To include population
cohort-share uncertainty, follow {cmd:jwselect_sa} with
{cmd:jwselect_sa_iw, postwindow(a b)} and then run {cmd:jwselect}.
That source forms one population IW target and the full joint covariance
with its cohort-specific preperiod restrictions, including the influence
of the estimated shares. Selection uses scalar target-level MD, rather than
averaging separately selected cohort estimates. See {help jwselect_sa_iw}.
Importing raw, unverified
{cmd:eventstudyinteract} or {cmd:reghdfe} results is not supported.

{title:Optional CS-DR adapter}

{p 4 4 2}
Install {cmd:csdid} separately using {browse "https://psantanna.com/csdid/":the authors' installation guide}.
The adapter targets the 2.0.1 API and requires version 2.0.1 or later in the 2.x
line; older SSC 1.x releases are rejected. The tested 2.0.1 source commit is
{cmd:f98bcea0691b1670c409bb7e3cc81d60705b4f84}. Record {cmd:csdid version} for reproducibility. No dependency is installed by
{cmd:jwselect}; the matrix interface needs no upstream estimator.

{p 4 4 2}
Run {cmd:csdid} with {cmd:method(dr) nevertreated base_period(universal) analytical pointwise pscoretrim(1)}.
The supported sample is an unweighted balanced panel with positive integer,
unit-spaced calendar time, actual never-treated controls coded zero, and plain
numeric covariates constant within unit. Omit extra clustering or use
{cmd:cluster(ivar)}; higher-level clusters are unsupported. Keep the original
estimation data. Call {cmd:jwselect} on unaggregated ATT(g,t) results, before
posting an aggregation. Bootstrap covariance, varying bases, not-yet-treated
controls, and unsupported samples are rejected.

{p 4 4 2}
{cmd:pscoretrim(1)} is required. The upstream default {cmd:pscoretrim(.995)}
is rejected even if no sample observations happen to be trimmed. Active
trimming changes the untrimmed target and can invalidate the DR union used
here. A retained or clipped propensity score is not a substitute for the
overlap assumptions. The composition-changing repeated-cross-section route
is supplied by {cmd:jwselect_rcs}, described below.

{p 4 4 2}
{cmd:vce(supplied)} is the only CS-DR covariance convention and is the default
for that source. The adapter retains the full analytical influence-function
covariance, including cross-cohort and cross-period terms. It reconstructs
structural-zero reference cells from the upstream cell metadata.
{cmd:baseline(cohort)} is required; an omitted {cmd:anticipation()} is inferred
from {cmd:csdid}'s value, and an explicit value must match it.

{title:Doubly robust repeated cross sections}

{p 4 4 2}
Run {cmd:jwselect_rcs y x1 x2, tvar(time) gvar(first_treat) cohort(g) target(s) baseline(b)}
before {cmd:jwselect}. The target is the ATT for cohort g's treated population
observed at date s. Cohort-date outcome means are standardized to that target's
covariate distribution; these standardized means are called transported means.
The distributions of cohort and covariates may change over time.
Each source cell uses an outcome regression and a
normalized target-to-source density ratio; the DR union requires either its
outcome regression or its ratio model to be correct. All target, nuisance
and normalization contributions enter the joint influence covariance.

{p 4 4 2}
{cmd:sampling(pooled)} (default) treats observations as an independent pooled
mixture. {cmd:sampling(stratified)} allows independent date-specific samples
with fixed date sample sizes and centers the complete influence within each
date. Both require representative independent samples, covariates unaffected
by treatment, overlap and the maintained target-standardized average
parallel-trends equality. Changing composition does not license arbitrary
outcome-selective sampling. Overlapping panels, survey weights, trimming and
extra clustering are unsupported. See {help jwselect_rcs} for source syntax,
optional preperiods and detailed assumptions.

{p 4 4 2}
Changing composition makes the moments target-specific. Refit the source
for a different target; separately fitted targets do not supply their
cross-target covariance. The scalar source supports manual, WBIC and
aggregate-focused selection, with the source's full sampling covariance.

{title:Optional nonlinear jwdid adapter}

{p 4 4 2}
{cmd:method(logit)} and {cmd:method(poisson)} are supported under the balanced,
unweighted, never-treated, unit-clustered and time-invariant-covariate contract
above. Other nonlinear methods, including {cmd:ppmlhdfe}, are not supported.
Omitting {cmd:restrictions()} or specifying {cmd:restrictions(index)} uses
the full conditional-index formulation. Under its correctly specified finite-dimensional conditional model
and identifying restrictions, it selects the entire lead coefficient block:
the intercept and every covariate interaction at a cohort/date. The ATT
remains on the outcome scale. Unit scores provide the full joint covariance
between ATT and index coefficients. {cmd:vce(unconditional)} includes
empirical target-covariate-distribution uncertainty and its cross terms;
{cmd:vce(conditional)} conditions on the observed covariates. Index
coefficients are reported in the original covariate basis.
Centered interactions and upstream {cmd:xasis} are supported.

{p 4 4 2}
In index mode each complete date block is one selection unit. With q covariates it has q+1
identified coefficients. Wald-BIC penalizes the restriction rank, not the number
of selected dates. Generated candidate families enumerate dates; {cmd:prematrix()}
selects complete blocks. Custom candidate masks have {cmd:e(K)} coefficient
columns in {cmd:e(hmap)} order and must assign the same bit to every component
of a date. Split, unidentified, or singular blocks are rejected.

{p 4 4 2}
Explicit {cmd:restrictions(placebo)} retains the legacy outcome-scale option:
joint {cmd:estat attgt} margins supply the ATT and average fitted-placebo
moments. Each optional date contributes one scalar restriction. This is a
different moment choice, not an equivalent implementation of the
conditional-index restrictions. A zero average placebo can hide offsetting
nonzero index coefficients. Identification and validity of these average
moments require their own substantive justification.
{cmd:restrictions()} is accepted only after supported nonlinear logit or Poisson fits.

{p 4 4 2}
A weaker target-specific average formulation is a possible additional
extension. It requires a separately established identifying map and joint
estimation/covariance argument. It is not supplied by the legacy placebo
option, and there is no {cmd:restrictions(average)} option.

{p 4 4 2}
Retain {cmd:baseline(cohort)} and the upstream normalization. If
{cmd:anticipation()} is omitted, the adapter uses the recorded {cmd:jwdid}
anticipation value minus one. An explicit incompatible value, or
{cmd:baseline(common)}, is rejected. The original linear/matrix rebasing
behavior is unchanged.

{title:Moment interpretation for CS-DR and nonlinear inputs}

{p 4 4 2}
The maintained estimate d is an outcome-scale ATT under the source's
identifying assumptions; optional h values are average fitted-placebo
moments, DR placebo contrasts, or full conditional-index lead blocks
under {cmd:restrictions(index)}. For a fixed candidate, the adjustment
is d - Cov(d,h) inv(Var(h)) h, using the same joint covariance throughout.
This is a Gaussian covariance-weighted moment adjustment. Its Wald--J identity
belongs to the canonical quadratic MD problem with unrestricted post effects
and h=0. It does not assert equality with a constrained nonlinear refit or
CS-DR re-estimation at alternative base periods.

{p 4 4 2}
With covariate interactions, zero average nonlinear fitted-placebo moments
can hold while index-pretrend coefficients are nonzero. They do not establish
conditional-index parallel trends or the proposed target-specific extension.
The index route
uses the full coefficient block for target-level MD; the separate
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
sampling-covariance scale. This equals the estimated asymptotic mean squared
error (AMSE) criterion divided by N; it gives the same minimizer. The squared bias correction is
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
After {cmd:jwselect_rcs} or {cmd:jwselect_sa_iw}, the source supplies one scalar
population target. Only that estimate is reported. Focused selection defaults
to {cmd:focus(aggregate)}; explicit {cmd:focus(cohort)} or {cmd:focus(both)} is
rejected. The restriction maps still use their actual cohort/date labels,
so manual masks, WBIC scopes and joint candidate families retain their meanings.

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

{p 4 4 2}
The November 16, 2026 manuscript's empirical workflow uses cohort-wise
all-subsets selection. Specify {cmd:wbicscope(cohort) candidates(all) maxcand(16384)}
to allow up to 14 eligible optional dates per cohort (2^14 candidates,
including the empty set). This differs from the computational defaults,
{cmd:candidates(nested) maxcand(4096)}. A larger family requires a larger
explicit cap; this choice does not change the final joint-MD step.

{title:Matched baseline, own-cohort, and joint MD comparison}

{p 4 4 2}
{cmd:estat compare} holds the selected restrictions, target definitions,
aggregation weights, and supplied covariance fixed. The baseline estimate
uses no optional restrictions. Own-cohort MD adjusts each cohort's effects
using only its own part of the selected mask. Joint MD adjusts those effects
using the full mask and cross-cohort covariance. No selection is rerun.
The aggregate variance for own-cohort MD retains cross-cohort covariance;
it is not the sum of isolated cohort variances.

{p 4 4 2}
By default the table follows {cmd:e(output_map)}. {cmd:output(#)} selects
one row of that map. {cmd:cells} instead compares the primary cohort-time
cell system and returns {cmd:r(postmap)}. The two options are mutually
exclusive. The scalar population sources {cmd:jwselect_rcs} and
{cmd:jwselect_sa_iw} do not supply a cohort decomposition and are rejected.

{p 4 4 2}
For manual and WBIC selection, all outputs use the same selected union.
For cohort-focused primary results, the joint comparison uses the union
of cohort selections; the own-cohort comparison uses each cohort's slice.
A separately aggregate-focused row uses its own selected joint mask for
both comparisons. Thus focused rows may use different masks. Comparing
these estimates does not extend focused inference guarantees to another
reporting target or to the new joint counterpart.

{p 4 4 2}
The columns of {cmd:r(table)} are {cmd:baseline se_baseline own se_own joint se_joint}
{cmd:own_minus_base joint_minus_base joint_minus_own}. All standard errors
are calculated by treating the selected restrictions as fixed. They are
descriptive diagnostics, not selection-adjusted inference or tests of
estimator differences.
The command does not produce confidence intervals or p-values. Existing
{cmd:e(ci)} and {cmd:e(pvalue)} remain attached to the original estimator
and inference method; do not transfer them to another comparison column.
Neither the variance comparison nor the direction of an estimate's change
guarantees greater accuracy in finite samples.

{p 4 4 2}
Saved matrices include {cmd:r(V_baseline)}, {cmd:r(V_own)}, {cmd:r(V_joint)},
{cmd:r(loading_baseline)}, {cmd:r(loading_own)}, {cmd:r(loading_joint)},
and {cmd:r(selection)} when K>0. The mask columns follow {cmd:e(hmap)}.
The map is {cmd:r(output_map)} for reported outputs and {cmd:r(postmap)}
for cells. {cmd:r(output)} records the requested output number (0 for all),
{cmd:r(cells)} is 1 for cell output, and {cmd:r(inference)} is
{cmd:fixed_selection_diagnostics}. The command preserves the active estimation results, data,
and RNG state. Copy the required r() matrices before running another
postestimation command.

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
The command treats supplied {cmd:aggweights()} as fixed. For population
Sun-Abraham shares, use {cmd:jwselect_sa_iw}, which adds their complete
influence contribution before selection. Do not substitute estimated shares
into {cmd:aggweights()} and interpret its covariance as population-IW inference. Signed
contrasts are outside this initial averaging interface.

{title:Confidence intervals}

{p 4 4 2}
Manual selection reports normal intervals under the fixed restrictions.
By default, Wald-BIC reports an ordinary pointwise normal interval treating the
selected restrictions as fixed. To request a marginal selection-adjusted interval, specify
{cmd:inference(wbic)} and explicitly choose one of the following methods.
Both retain the ordinary WBIC point estimate. Use {cmd:estat inference}
to inspect the intervals and method-specific diagnostics.

{p 4 4 2}
{cmd:wbicmethod(bootstrap)} implements centered Gaussian bootstrap reselection:
each draw from the fitted Gaussian experiment repeats preperiod selection and
MD estimation, using a centered, studentized error statistic. This numerically approximates
a pointwise asymptotic construction; no higher-order accuracy improvement
is asserted. A fixed number of draws leaves quantile simulation error and
does not by itself satisfy a vanishing numerical-error condition. Increase
draws and examine numerical sensitivity. This is not a uniform or finite-sample
coverage guarantee.
{cmd:wbicengine(integrated)} is the default; {cmd:wbicengine(direct)} uses direct
Gaussian simulation. {cmd:reps()} and {cmd:seed()} set the simulation budget and seed.

{p 4 4 2}
{cmd:wbicmethod(mc)} uses the smallest candidate probability level that passes the
McCloskey rejection-probability bound in the fitted Gaussian experiment.
{cmd:mcsteps(9)} sets the number of candidate levels. {cmd:mcrefine} optionally
searches for a lower level satisfying this bound. The calibration allowance is
{cmd:delta=mcerr()/2} with positive rejection-probability target
{cmd:1-level()/100-delta}. McCloskey does not use {cmd:gamma()}.
When no search level passes, the default {cmd:mcfallback(terminal)} returns
the analytic terminal interval in the same family.
Its finite half-width is evaluated on the normal-score scale and contains the
baseline interval at the adjusted error level. {cmd:mcfallback(error)} stops
when the search fails. Status 7 identifies a finite analytic terminal. Status 8
denotes an exceptional whole-line return, displayed as [-inf, +inf]. Because
Stata matrices have no signed infinity representation, its {cmd:e(ci)} endpoints
are missing and {cmd:e(mc_terminal_diagnostics)} records {cmd:unbounded=1}.

{p 4 4 2}
The uniform asymptotic coverage result requires the stated regularity conditions
and rejection-probability bounds that remain conservative after numerical
approximation and rounding. These bounds are evaluated in double precision
with numerical padding. This implementation does not supply a formal
interval-arithmetic certificate for all floating-point operations.
The retained {cmd:mcmode(certified)} option and diagnostic labels refer to
the implemented Gaussian calibration route. They are compatibility labels,
not a claim that conservative bounds have been formally verified for every
floating-point operation. Monte Carlo error allocation and numerical roundoff control
are distinct requirements.

{p 4 4 2}
Read the reported endpoints from {cmd:e(ci)} or {cmd:e(table)}. The adjusted
methods report no adjusted p-values. Displayed standard errors and {cmd:e(V)}
are calculated by treating the selected restrictions as fixed.
{cmd:inferoutputs()} limits adjusted inference
to specified rows of {cmd:e(output_map)}.

{title:Focused inference}

{p 4 4 2}
For focused selection, {cmd:inference(auto)} uses
{cmd:inference(reselect)}. At every trial nuisance value, common joint Gaussian
draws rerun the actual focused candidate rule. The selected estimation error
includes the nuisance mean. Selection and reporting targets are distinguished:
aggregate-focused cohort rows use the aggregate's selection rule, and custom
weighted aggregates of cohort-focused estimates retain each cohort's original
selection rule. Aggregates are formed within each draw before taking quantiles,
preserving cross-cohort covariance.

{p 4 4 2}
When selection and reporting targets differ, the cross-target extension
requires additional regularity, including the stated condition on limiting
ties of selection scores and reporting adjustments. The software does not
automatically verify that condition. The existence of a computed interval
does not establish it. See the full guide's cross-target extension before
interpreting such intervals as an application of its coverage result.

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
All displayed standard errors are labeled {cmd:se_fixed}. Under focused selection
they are diagnostics calculated by treating the selected restrictions as fixed,
not standard errors for a normal approximation to the selected estimator.
Accordingly, focused results post e(b) without e(V).
Use e(ci) and e(pvalue). e(V_fixed) and e(V_cells) preserve the joint covariance
calculated by treating the selected restrictions as fixed for diagnostics;
do not pass them to conventional
{cmd:lincom}, {cmd:test}, or {cmd:margins} as selection-adjusted covariance.
No single covariance matrix characterizes the focused estimator's nonnormal
selection distribution. Observation-level {cmd:predict} is deliberately undefined.

{title:Requested outputs, caching, and prepared jobs}

{p 4 4 2}
{cmd:inferoutputs(numlist)} identifies output ROW NUMBERS in e(output_map),
not cohort dates or event times. It supports focused selection with
{cmd:inference(none)}, {cmd:auto}, or {cmd:reselect}, and WBIC selection with
{cmd:inference(wbic)}. Without it, inference is
requested for every output. Unrequested outputs retain point estimates and
covariance calculated by treating the selected restrictions as fixed, but have
missing intervals and p-values. Use this
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
{cmd:e(scalar_target)} identifies a scalar population source and
{cmd:e(target_label)} describes its target. For scalar sources, the sole
{cmd:e(postmap)} row is a representative label, not a disaggregated cell;
{cmd:estat cells} is unavailable. The joint canonical vector contains that
one target followed by all optional moments. Its complete covariance includes
the source-specific composition or population-share contribution.

{p 4 4 2}
Key matrices: e(b), e(V) for manual/WBIC only, e(V_fixed), e(table), e(ci),
e(pvalue), e(output_map), e(baselines), e(postmap), e(hmap), e(selection),
e(cell_selection), e(b_cells), e(V_cells), e(weights), e(scores),
e(candidates), e(joint_candidates), e(z), e(V_z), e(loading), e(targets).
Numerical reselection also stores e(inference_diagnostics) and
e(inference_extrema).
The centered bootstrap stores {cmd:e(wbicengine)} and {cmd:e(wbic_pointwise_diagnostics)}.
McCloskey stores {cmd:e(mc_diagnostics)}, {cmd:e(mc_profile_diagnostics)},
{cmd:e(mc_terminal_diagnostics)} and, when calibration is attempted,
{cmd:e(mc_calibration)}. Terminal diagnostic columns are {cmd:output}, {cmd:used},
{cmd:unbounded}, {cmd:s0}, {cmd:k_safe}, {cmd:half_width}, {cmd:epsilon} and
{cmd:baseline_halfwidth}. The terminal uses {cmd:epsilon=1e-6}. Its {cmd:pstar}
is stored missing because its probability may round to one. Retain {cmd:k_safe}
when interpreting that return.
The stored diagnostic column names {cmd:risk_upper} and {cmd:risk_target}
remain unchanged for compatibility; they denote the rejection-probability
upper bound and its calibration target.
The {cmd:mc_diagnostics} status is 0 for calibration, 1 for an unrequested output,
2 for no-restriction normal inference, 7 for the analytic terminal and 8 for
the exceptional whole-line return. The {cmd:rule_fallback} column identifies
the observation-domain guard within the McCloskey rule.
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
The source estimator is matrix, csdid, jwdid, eventstudyinteract, jwselect_rcs,
or jwselect_sa_iw. e(source_version) is empty
when the source does not expose a verified version, including jwdid; retain
{cmd:which jwdid} output separately. e(version) is {cmd:0.9.1} and
e(build) is {cmd:0.9.1-20261116}. The exact build identifies the numerical
modules and guards against mixing releases in one session.
Start a fresh Stata session and recreate prepared jobs after upgrading.

{p 4 4 2}
{cmd:e(wbicmethod)} records the adjusted method. For McCloskey, inspect
the calibration, terminal and guard flags by output with {cmd:estat inference}.

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

{phang2}{cmd:csdid y x1 x2, ivar(id) time(year) gvar(first_treat) method(dr) nevertreated base_period(universal) analytical pointwise pscoretrim(1) cluster(id)}
{phang2}{cmd:estimates store CS_SOURCE}
{phang2}{cmd:jwselect, select(wbic) postwindow(0 2)}
{phang2}{cmd:estat compare}
{phang2}{cmd:matrix comparison = r(table)}
{phang2}{cmd:estat compare, cells}
{phang2}{cmd:estimates restore CS_SOURCE}
{phang2}{cmd:jwselect, select(manual) postwindow(0 2)}

{phang2}{cmd:jwdid y x1 x2, ivar(id) tvar(year) gvar(first_treat) never method(logit) cluster(id)}
{phang2}{cmd:estimates store NL_SOURCE}
{phang2}{cmd:jwselect, select(wbic) restrictions(index)}
{phang2}{cmd:estimates restore NL_SOURCE}
{phang2}{cmd:jwselect, select(manual)}

{p 4 4 2}
For an exponential conditional mean use {cmd:method(poisson)} in the source fit.
The explicit {cmd:restrictions(index)} above is also the nonlinear default.
Restore the stored source fit before each
new direct specification. The following example retains the linear workflow.


{phang2}{cmd:jwdid y x1 x2, ivar(id) tvar(year) gvar(first_treat) never method(regress)}
{phang2}{cmd:estimates store original}
{phang2}{cmd:jwselect, select(wbic) wbicscope(joint) candidates(nested)}

{p 4 4 2}
Restore the original source fit before each comparison.
{phang2}{cmd:estimates restore original}
{phang2}{cmd:jwselect, inference(wbic) wbicmethod(bootstrap) reps(19999) seed(12345)}
{phang2}{cmd:estat inference}
{phang2}{cmd:estimates restore original}
{phang2}{cmd:jwselect, inference(wbic) wbicmethod(mc)}
{phang2}{cmd:estat inference}
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
A self-contained matrix-input example is in {cmd:examples/matrix_demo.do}
in the downloaded package.

{p 4 4 2}
For an explicit all-subsets, cohort-wise WBIC comparison with centered
Gaussian reselection, retain the same prespecified event window in all
steps. This template illustrates the manuscript's workflow; its window,
baseline, anticipation and covariance must be chosen for the application.
It does not automatically reproduce the manuscript's empirical targets.

{phang2}{cmd:estimates restore original}
{phang2}{cmd:jwselect, select(wbic) wbicscope(cohort) candidates(all) maxcand(16384) postwindow(0 2) inference(wbic) wbicmethod(bootstrap) reps(49999) seed(12345)}
{phang2}{cmd:estat compare}
{phang2}{cmd:matrix comparison = r(table)}
{phang2}{cmd:estat inference}
