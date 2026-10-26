*! jwselect 0.5.0: nonlinear DiD with optional jwdid dependency.
*! Run from the extracted package directory. Replaces data and estimates.
*! Installs nothing; jwdid and its dependencies must already be available.
version 16.1
set more off
adopath ++ "."
capture which jwdid
if _rc {
    display as error "This example requires jwdid and its dependencies; see README.md."
    exit 111
}

clear all
set seed 681907
set obs 4800
generate long id = ceil(_n/8)
generate byte year = mod(_n-1,8)+1
generate byte first_treat = cond(mod(id,3)==0,0,cond(mod(id,3)==1,4,6))
generate double x1 = sin(id/9)
generate double x2 = cos(id/13)
generate byte treated = first_treat>0 & year>=first_treat
generate double index = -.6 + .08*year + .25*x1 - .15*x2 + .4*treated
generate byte ybinary = runiform() < invlogit(index)
generate long ycount = rpoisson(exp(index))

* Logit: jwselect consumes outcome-scale margins, not logit coefficients.
jwdid ybinary x1 x2, ivar(id) tvar(year) gvar(first_treat) ///
    never method(logit) cluster(id)
estimates store LOGIT_SOURCE
jwselect, select(manual) postwindow(0 2)
estimates store LOGIT_BASE
estimates restore LOGIT_SOURCE
jwselect, select(wbic) candidates(nested) postwindow(0 2)
estat selection
estat pretrend
estimates store LOGIT_WBIC

* Poisson uses an exponential conditional mean under the same adapter contract.
jwdid ycount x1 x2, ivar(id) tvar(year) gvar(first_treat) ///
    never method(poisson) cluster(id)
estimates store POISSON_SOURCE
jwselect, select(manual) postwindow(0 2)
estimates store POISSON_BASE
estimates restore POISSON_SOURCE
jwselect, select(wbic) candidates(nested) postwindow(0 2)
estat selection
estimates store POISSON_WBIC

* vce(unconditional) is the default. The upstream cohort baseline is retained;
* compatible anticipation is inferred when omitted. Do not request rebasing.
* Optional restrictions are zero average outcome-scale placebo effects.
* Gaussian MD adjustment is not a constrained nonlinear refit and may leave
* natural outcome bounds; jwselect does not clip the results.
ereturn list
