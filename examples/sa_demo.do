*! jwselect 0.7.0: Sun-Abraham, baseline and selected preperiods.
*! Synthetic data; run from the extracted package root. Clears data/estimates.
version 16.1
clear all
set more off
adopath ++ "."
foreach dependency in eventstudyinteract avar reghdfe ftools {
    capture which `dependency'
    if _rc {
        display as error "Install optional dependency `dependency' first; see README.md."
        exit 499
    }
}
set seed 29102026
set obs 600
generate long id = _n
generate byte first_treat = cond(mod(id,3)==0,0,cond(mod(id,3)==1,4,6))
generate double alpha = rnormal()
expand 8
bysort id: generate byte time = _n
generate double y = alpha + .2*time + rnormal() + ///
    (first_treat>0 & time>=first_treat)*(1+.15*(time-first_treat))

jwselect_sa y, ivar(id) tvar(time) gvar(first_treat)
estimates store SA_SOURCE

* Baseline only; equal cohort weights over event times 0,1,2.
jwselect, select(manual) postwindow(0 2)
estimates store SA_BASE

* Joint Wald-BIC retains effect/placebo cross-covariances.
estimates restore SA_SOURCE
jwselect, select(wbic) wbicscope(joint) candidates(nested) postwindow(0 2)
estat selection
estat pretrend
estimates store SA_WBIC

* Fixed, prespecified extra restrictions (calendar dates).
matrix PRE = (4,2 \ 6,4)
estimates restore SA_SOURCE
jwselect, select(manual) prematrix(PRE) postwindow(0 2)

* Focused selection only. This does NOT produce selection-adjusted CIs.
estimates restore SA_SOURCE
jwselect, select(focused) focus(both) postwindow(0 2) inference(none)
estat scores

* For selection-adjusted CIs, restore SA_SOURCE and use inference(reselect).
* The default simulation/search settings are much more costly; inspect
* estat inference and numerical stability before reporting intervals.
