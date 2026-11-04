{smcl}
{* *! jwselect 0.9.0 04nov2026}{...}
{title:jwselect_sa_iw -- Population-share Sun-Abraham target}

{title:Syntax}
{p 8 12 2}
{cmd:jwselect_sa_iw}, {opt postwindow(a b)} [{opt cohorts(numlist)}
{opt baseline(cohort|common)} {opt anticipation(#)} {opt level(#)}]

{title:Workflow}
{p 8 12 2}{cmd:jwselect_sa y, ivar(id) tvar(year) gvar(first_treat)}
{p 8 12 2}{cmd:estimates store SA_SOURCE}
{p 8 12 2}{cmd:jwselect_sa_iw, postwindow(0 2)}
{p 8 12 2}{cmd:estimates store IW_SOURCE}
{p 8 12 2}{cmd:jwselect, select(wbic) candidates(all)}

{p 4 4 2}
This command converts a validated, signed {help jwselect_sa} fit into one scalar
population target with its full joint covariance with cohort-specific preperiod
restrictions. It averages event times a,...,b equally within each supported
cohort, then averages cohorts with estimated population shares among included
cohorts. Use {cmd:postwindow(e e)} for the usual population IW effect at one
exposure time. Optional moments retain every individual cohort/date contrast;
they are never replaced by an averaged zero lead.

{p 4 4 2}
By default the cohort set consists of cohorts supporting the entire requested
window. {cmd:cohorts()} fixes a subset explicitly; an unsupported requested
cohort is an error. The supported cohort set is held fixed in the theory.
This does not provide inference for data-dependent discovery of cohort support.
The source requires the balanced, unweighted, no-covariate, never-treated design
of {cmd:jwselect_sa}; no new regression is fitted.

{title:Target and covariance}
{p 4 4 2}
Let d_g be the baseline-identified average effect in the requested window and
w_g=P(G=g)/P(G in K). The target is sum_g w_g d_g. Its influence is the weighted
sum of cell-contrast influences plus
1{G in K}(d_G-target)/P(G in K). The latter term accounts for estimating the
population cohort shares. The complete joint panel influence covariance is
computed directly, retaining all cross terms. It has no regression
degrees-of-freedom multiplier. Do not add a further share correction.

{p 4 4 2}
Baseline validity, no anticipation at included dates, independent panels,
positive cohort/control probabilities, and the relevant moment and covariance
regularity conditions remain required. Numerical selection-adjusted coverage
retains the numerical qualifications in {help jwselect}. This is target-level
MD for the scalar population target; it is not the same operation as averaging
separately selected cohort effects with estimated shares.

{p 4 4 2}
The existing direct {cmd:jwselect_sa} route retains fixed-weight reporting.
Passing estimated shares to {cmd:aggweights()} in that route does not substitute
for this command. Restore {cmd:SA_SOURCE} to build another population target;
restore {cmd:IW_SOURCE} before comparing selection methods for the same target.

{title:Saved results}
{p 4 4 2}
{cmd:e(jw_d)} is the scalar baseline estimate; {cmd:e(jw_h)} the row of preperiod
contrasts; {cmd:e(jw_V)} is the sampling covariance ordered (d,H).
{cmd:e(jw_hmap)} gives cohort and calendar date of each restriction;
{cmd:e(jw_bases)} gives included cohorts and maintained baselines.
Empty H matrices are omitted when no optional restrictions exist.
{cmd:e(jw_nunits)} is the independent panel count.

{p 4 4 2}
{cmd:e(iw_cohorts)} has columns cohort, N, population share estimate, normalized
aggregation weight, and baseline cohort-window effect.
{cmd:e(iw_targetweights)} lists each target cohort/calendar cell and its weight.
{cmd:e(iw_share_variance)} and {cmd:e(iw_fixed_variance)} report the share and
fixed-weight pieces of the baseline scalar variance. {cmd:e(b)} and {cmd:e(V)}
contain the scalar estimate and its full variance. The signed data and e(sample)
are retained for validation by {cmd:jwselect}.
