{smcl}
{* *! jwselect 0.9.0 04nov2026}{...}
{title:Title}
{p 4 4 2}{cmd:jwselect_rcs} -- Doubly robust DiD inputs for repeated cross sections allowing changing composition

{title:Syntax}
{p 8 12 2}
{cmd:jwselect_rcs} {it:outcome} [{it:covariates}] {ifin},
{cmd:tvar(}{it:time}{cmd:)} {cmd:gvar(}{it:first_treatment}{cmd:)}
{cmd:cohort(}{it:g}{cmd:)} {cmd:target(}{it:s}{cmd:)}
{cmd:baseline(}{it:b}{cmd:)}
[{cmd:preperiods(}{it:numlist}{cmd:)} {cmd:gap(}{it:integer}{cmd:)}
{cmd:psvars(}{it:varlist}{cmd:)} {cmd:orconstant} {cmd:psconstant}
{cmd:sampling(pooled|stratified)}
{cmd:iterate(}{it:integer}{cmd:)} {cmd:level(}{it:#}{cmd:)}]

{title:Target and identifying assumptions}
{p 4 4 2}
Each invocation fixes one cohort g and one post-treatment date s. Never-treated
observations have gvar=0. The target is the ATT for the treated population
observed at date s, conditional on G=g and R=s. All preperiod and control
conditional means are standardized to that same target covariate distribution.
The distribution of (G,X) can change over dates; stationarity is not imposed.

{p 4 4 2}
The identifying baseline requires the target's untreated mean at s to equal
the target-standardized sum of the treated baseline mean and the control change
from b to s. This is the average form of cross-population conditional parallel
trends. No anticipation at every included treated preperiod must hold after
standardization to the target covariate distribution. Conditional no
anticipation within each treated preperiod source cell is sufficient;
zero average anticipation under that source cell's own covariate distribution
alone need not remain zero after transport. This requirement applies to the
maintained baseline and each included optional preperiod. It does not permit
arbitrary outcome-selective sampling. Maintain representative independent
cross-sectional observations, unaffected covariates, and common support
of every source cell over the target population. Panel overlap, survey weights,
and dependence across observations are not supported by this first-stage command.

{p 4 4 2}
Changing composition generally makes restrictions target-specific. Refit for
another target date; do not combine separate runs as if they used one common
preperiod vector or infer their cross-target covariance from marginal outputs.

{title:Estimation and covariance}
{p 4 4 2}
For each source cell j=(G,R), the command fits a linear outcome regression in
that cell and a binary logit comparing the target cell A=(g,s) with j. The
exponentiated logit index estimates the density ratio P(A|X)/P(j|X). Its
transported mean equals the target average of the fitted outcome regression
plus the normalized source-weighted residual average. The target cell's own
mean is its observed average and requires no target outcome regression.

{p 4 4 2}
Each transported mean is consistent if its outcome regression is correct OR
its target-to-source density ratio is correct, under overlap and regularity.
This union condition must cover every cell used in the target and optional
contrasts. Correct source control changes alone do not suffice under unrestricted
composition changes. No semiparametric efficiency claim is made under arbitrary
parametric misspecification.

{p 4 4 2}
The covariance uses the complete observation-level influence, including outcome
and propensity fitting, normalized weights, empirical target-distribution
uncertainty, and shared-observation cross-covariances. It remains a full
parametric sandwich when only one nuisance component is correct. There is no
degrees-of-freedom multiplier. N is the number of eligible pooled observations,
not a panel count or a number of treated observations.

{p 4 4 2}
{cmd:sampling(pooled)}, the default, treats observations as independent draws
from a date-mixture population. {cmd:sampling(stratified)} instead conditions on
fixed date-specific sample sizes, with independent representative draws within
each date and positive limiting date fractions. It centers the COMPLETE fitted
influence within each date before forming cross-products. Estimates are unchanged;
only the covariance changes. This is date-stratified sampling, not support for
survey weights, clustered/overlapping panels, or arbitrary strata.

{p 4 4 2}
Trimming and clipping are not implemented or permitted. Singular fits,
nonconvergence, nonfinite density ratios, or fitted pairwise propensities
numerically within 1e-12 of zero or one stop with an error. This refusal is
not trimming: no modified estimator is silently returned. These empirical
checks do not establish population overlap or uniform regularity.

{title:Options}
{p 4 4 2}
{cmd:preperiods()} lists optional dates no later than b and on the specified
time grid. Including b is harmless; it is retained once as the structural zero.
With no optional dates, jwselect performs baseline-only inference.

{p 4 4 2}
{cmd:gap()} defaults to 1. Dates and cohort must be positive integers, b<g<=s,
and g-b and s-g must be multiples of gap. The implied anticipation allowance
is (g-b)/gap-1; all optional dates precede g by more than that allowance.

{p 4 4 2}
Covariates are ordinary numeric variables. Generate factor indicators and
transformations explicitly; the command supplies an intercept. By default both
nuisance models use the listed covariates. {cmd:psvars()} supplies a different
propensity-regression list. {cmd:orconstant} or {cmd:psconstant} restricts that
working model to an intercept; the corresponding other model must be correct
for the DR union guarantee. {cmd:iterate()} defaults to 100.

{title:Selection}
{p 4 4 2}
After fitting, invoke {cmd:jwselect} with manual, WBIC, or focused selection.
The scalar source fixes the population target and baseline. Its preperiod
restrictions are differences of target-standardized treated-control gaps.
The common joint covariance is used for target-level MD. The ordinary pointwise,
focused-local, and McCloskey uniform/numerical qualifications still apply.

{phang2}{cmd:jwselect_rcs y x1 x2, tvar(year) gvar(first_treat) cohort(2005) target(2008) baseline(2004) preperiods(2001 2002 2003)}
{phang2}{cmd:estimates store source}
{phang2}{cmd:jwselect, select(wbic) candidates(all)}
{phang2}{cmd:estimates restore source}
{phang2}{cmd:jwselect, select(focused) candidates(all) inference(none)}

{title:Saved results}
{p 4 4 2}
e(bgt), e(Vgt), e(cells) hold optional gap contrasts, the structural-zero
baseline, and the target ATT with joint sampling covariance. e(transported_means)
lists source cohort, date, and transported mean; e(transported_V) is their
joint sampling covariance. e(rcs_diagnostics) contains source count, pairwise
propensity extrema, maximum unnormalized density ratio, and score residual.
The direct target row has missing propensity diagnostics because it needs no
propensity or outcome fit. e(jw_d), e(jw_h), e(jw_V), and related scalar-schema
metadata supply the selection adapter; when e(jw_K)=0, h and hmap are absent.
e(rcs_sampling) records pooled or stratified sampling.

{title:Reference}
{p 4 4 2}
Sant'Anna and Xu, {browse "https://arxiv.org/abs/2304.13925":Difference-in-Differences with Compositional Changes}.
The implemented multivariate preperiod construction transports every cell to
one fixed post-treated target; its finite-dimensional nuisance sandwich is
derived in the accompanying manuscript. The code does not implement the paper's
nonparametric or cross-fitting estimators.
