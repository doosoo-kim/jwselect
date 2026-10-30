*! jwselect 0.7.0: full conditional-index blocks, logit and Poisson.
*! Run from the extracted package root. Replaces data and estimates.
*! No dependencies are installed or updated by this example.
version 16.1
clear all
set more off
adopath ++ "."
capture which jwdid
if _rc {
    display as error "Install optional jwdid and its dependencies first; see README.md."
    exit 111
}
set seed 703018
set obs 5400
generate long id=ceil(_n/6)
generate byte time=mod(_n-1,6)+1
generate byte group=cond(mod(id,3)==0,0,cond(mod(id,3)==1,3,5))
generate double x1=sin(id/7)
generate double x2=cos(id/11)
generate byte treated=group>0 & time>=group
generate double eta=-.5+.08*time+.2*x1-.15*x2 ///
    +.3*treated+.12*treated*x1+.4*(group==5 & time==1)*x1
generate byte binary=runiform()<invlogit(eta)
generate long count=rpoisson(exp(eta))

foreach link in logit poisson {
    local outcome binary
    if "`link'"=="poisson" local outcome count
    jwdid `outcome' x1 x2, ivar(id) tvar(time) gvar(group) ///
        never method(`link') cluster(id)
    estimates store NL_SOURCE

    * Baseline-only outcome-scale ATT; no optional restrictions selected.
    jwselect, restrictions(index) select(manual) postwindow(0 1)
    estimates store BASE_`link'

    * Two covariates: three coefficients in each eligible date block.
    estimates restore NL_SOURCE
    jwselect, restrictions(index) select(wbic) wbicscope(joint) ///
        candidates(all) postwindow(0 1)
    estat design
    estat selection
    estat pretrend
    estimates store WBIC_`link'

    * Focused point estimates. To compute reselection intervals, restore
    * NL_SOURCE and change inference(none) to inference(reselect).
    estimates restore NL_SOURCE
    jwselect, restrictions(index) select(focused) focus(both) ///
        candidates(nested) postwindow(0 1) inference(none)
    estimates store FOCUSED_`link'
}

* Default vce(unconditional) includes covariate-distribution uncertainty.
* WBIC intervals are pointwise/oracle, not uniformly valid after selection.
* The final estimator is outcome-scale target-level MD; it is not a
* constrained nonlinear refit and does not preserve natural ATT bounds.
