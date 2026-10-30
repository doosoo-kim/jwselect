{smcl}
{* *! jwselect 0.7.0 30oct2026}{...}
{title:jwselect_sa -- Fit validated Sun-Abraham inputs for jwselect}

{title:Syntax}
{p 8 12 2}
{cmd:jwselect_sa} {it:outcome} {ifin}, {opt ivar(id)} {opt tvar(time)} {opt gvar(first_treat)} [{opt level(#)}]

{title:Purpose}
{p 4 4 2}
This fitting wrapper calls the optional {cmd:eventstudyinteract} command and
retains its underlying cohort-event-time coefficients and FULL joint regression
covariance. It constructs every observed relative-time indicator except -1.
It signs the original sample so {help jwselect} can check that the data have
not changed. Raw {cmd:eventstudyinteract} or unrelated {cmd:reghdfe} results
cannot be imported automatically: their event map is not preserved reliably.

{title:Supported design}
{p 4 4 2}
One numeric outcome; balanced unweighted panel; independent units with clustering
at {cmd:ivar()}; positive integer unit-spaced calendar times; treatment cohort
constant within unit; never-treated controls coded {cmd:gvar=0}. Every treated
cohort must have an observed g-1 baseline and at least one treated observation.
Use the estimation {cmd:if}/{cmd:in} sample to define a balanced analysis panel.
Each treatment cohort and the never-treated group require at least two units.

{p 4 4 2}
This is the unbinned, no-covariate, never-treated specification proved in
draft v46_12, Appendix {it:Sun-Abraham} and its implementation subsection.
Covariates, weights, bins, extra omitted leads, last-treated controls, repeated
cross sections, and higher-level or multiway clusters are outside this wrapper.
The treatment is assumed binary and absorbing; {cmd:gvar()} records its first
treatment date. These substantive assumptions are not tested by the wrapper.

{title:Optional dependencies}
{p 4 4 2}
Install {cmd:eventstudyinteract}, {cmd:avar}, {cmd:reghdfe}, and {cmd:ftools}
separately. See {browse "https://github.com/lsun20/EventStudyInteract":the author's repository}.
The compatibility target is author source version 0.1, 24 January 2022.
The wrapper never installs or updates dependencies. Other jwselect workflows
do not require them. See the package README for installation commands.

{title:Workflow}
{p 8 12 2}{cmd:jwselect_sa y, ivar(id) tvar(year) gvar(first_treat)}
{p 8 12 2}{cmd:estimates store SA_SOURCE}
{p 8 12 2}{cmd:jwselect, select(wbic) candidates(nested) postwindow(0 2)}
{p 8 12 2}{cmd:estat selection}
{p 8 12 2}{cmd:estat pretrend}
{p 8 12 2}{cmd:estimates restore SA_SOURCE}
{p 8 12 2}{cmd:jwselect, select(manual) postwindow(0 2)}

{p 4 4 2}
Restore the source before each new jwselect specification. Manual, cohort or
joint Wald-BIC, focused selection, and prepared worker/collector inference
use the existing jwselect options. Use {cmd:prematrix()} to specify manual
restrictions. All supported post cells are retained before selecting the
reporting window. The estimation sample and source results are restored on
failure; temporary event indicators are removed and RNG state is preserved.

{title:Normalization and covariance}
{p 4 4 2}
The fitted source omits event -1 only. jwselect uses linear contrasts
beta[g,t-g]-beta[g,b-g] for its maintained baseline b. Thus
{cmd:baseline(common)} and {cmd:anticipation()} can choose another supported
baseline without refitting the source. Their identifying assumptions remain
the user's responsibility. Only the omitted g-1 source cell is reconstructed
as a structural zero; unsupported interactions are never treated as estimated
zeros.

{p 4 4 2}
Use {cmd:vce(supplied)}, the default after {cmd:jwselect_sa}. The adapter uses
the entire unit-clustered regression covariance, including cross-cohort and
cross-period entries and its finite-sample adjustment. It does not use
{cmd:e(V_interact)}, which contains individual variances only, or the
aggregated {cmd:e(V_iw)}.

{title:Aggregation scope}
{p 4 4 2}
The reported jwselect targets use fixed equal-cohort, equal-cell, or custom
cell weights. They are not automatically the native population IW target.
Using estimated cohort shares in {cmd:aggweights()} treats those numbers as
fixed and does NOT add their sampling uncertainty. Draft v46_12 equation
(sa_if) supplies an additional share-influence term for population IW;
that target is not implemented in this release. Native IW matrices retained
from the source are diagnostics, not inputs to jwselect's selector.

{title:Saved results}
{p 4 4 2}
{cmd:level()} sets the source-cell confidence level. A bare {cmd:jwselect_sa}
replays its summary; replay accepts no options. {cmd:predict} is disabled.

{p 4 4 2}
{cmd:e(cmd)} is {cmd:jwselect_sa}. {cmd:e(b)} and {cmd:e(V)} retain the full
underlying source regression; temporary coefficient names are not a public
prediction interface. {cmd:e(sa_bgt)}, {cmd:e(sa_Vgt)}, and {cmd:e(sa_cells)}
hold the normalized cohort-time vector, full covariance, and cohort/calendar
map consumed by the adapter. The original source's {cmd:e(b_interact)},
{cmd:e(V_interact)}, {cmd:e(b_iw)}, {cmd:e(V_iw)}, and {cmd:e(ff_w)} remain
available. See {help jwselect} for results after selection.

{p 4 4 2}
{cmd:e(sa_Vif)} is the plug-in iid-panel influence covariance before the
regression finite-sample factor {cmd:e(sa_vce_factor)}. It is a validation
diagnostic; selection uses {cmd:e(sa_Vgt)}. The wrapper verifies direct
cell-mean equivalence and joint covariance equivalence up to that factor,
and rejects singular or numerically degenerate identified covariance.

{p 4 4 2}
The wrapper is for fitting and saving selection inputs. Do not use its raw
temporary-regressor coefficients as aggregate ATTs or use unrelated regression
postestimation on them. Subsequent inference is through {cmd:jwselect} and
its documented {cmd:estat} subcommands.
