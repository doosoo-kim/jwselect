*! jwselect 0.5.0: CS-DR with optional csdid 2.0.1+ dependency.
*! Run from the extracted package directory. Replaces data and estimates.
*! Installs nothing; see README.md for the pinned csdid installation.
version 16.1
set more off
adopath ++ "."
capture which csdid
if _rc {
    display as error "This example requires csdid 2.x, version 2.0.1 or later; see README.md."
    exit 111
}

clear all
set seed 481702
set obs 4800
generate long id = ceil(_n/8)
generate byte year = mod(_n-1,8)+1
generate byte first_treat = cond(mod(id,3)==0,0,cond(mod(id,3)==1,4,6))
generate double x1 = sin(id/9)
generate double x2 = cos(id/13)
generate double y = .02*id + .15*year + .03*x1*year - .2*x2 + ///
    (.7 + .3*x1)*(first_treat>0 & year>=first_treat) + .8*rnormal()

* Explicit options retain a common upstream baseline and analytical covariance.
csdid y x1 x2, ivar(id) time(year) gvar(first_treat) ///
    method(dr) nevertreated base_period(universal) ///
    analytical pointwise cluster(id)
estimates store CS_SOURCE

* Baseline only. vce(supplied) is the automatic csdid-input convention.
jwselect, select(manual) postwindow(0 2)
estimates store CS_BASE

* Restore the original source before every new direct specification.
estimates restore CS_SOURCE
jwselect, select(wbic) candidates(nested) postwindow(0 2)
estat selection
estat pretrend
estimates store CS_WBIC

estimates restore CS_SOURCE
jwselect, select(wbic) wbicscope(joint) candidates(nested) postwindow(0 2)
estat selection
estimates store CS_JOINT

* Fixed equal-cohort aggregation need not match csdid's default aggregation.
* Selected restrictions concern mean outcome-scale placebo effects.
* This adjustment does not re-estimate csdid with alternative base periods.
ereturn list
