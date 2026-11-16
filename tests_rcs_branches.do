version 17
clear all
set more off
capture log close _all
log using "./tests_rcs_branches.log",text replace
adopath ++ "."
set seed 62712026
set obs 24000
gen double x=rnormal()
gen double u=runiform()
gen double total=0
forvalues j=1/8 {
    gen double pr`j'=exp(.08*`j'+.12*(`j'-4.5)*x)
    replace total=total+pr`j'
}
gen double cum=0
gen int cell=.
forvalues j=1/8 {
    replace cum=cum+pr`j'/total
    replace cell=`j' if missing(cell) & u<cum
}
gen int time=mod(cell-1,4)+1
gen int cohort=cond(cell<=4,4,0)
gen double y=1+.5*x+.3*time+.8*(cohort==4)+2*(cohort==4 & time==4)+rnormal()
quietly jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) preperiods(1 2) psconstant
assert abs(_b[g4_t4]-2)<.15
matrix ORB=e(b)
matrix ORV=e(V)
do ./tests_rcs_reference.mata
// Outcome correct, propensity deliberately misspecified; estimate remains centered.
di "RCS_OR_CORRECT_PS_WRONG_OK " _b[g4_t4]
replace y=y+.6*x^2
quietly jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) preperiods(1 2)
assert abs(_b[g4_t4]-2)<.15
// Outcome model now misses x^2; true pairwise log odds remain linear in x.
mata: rcs_test_reference()
assert scalar(rcs_reference_error)<1e-8
assert scalar(rcs_reference_b_error)<1e-8
di "RCS_PS_CORRECT_OR_WRONG_OK " _b[g4_t4]
di "RCS_UNION_SANDWICH_REFERENCE_OK " scalar(rcs_reference_error)
matrix POOLED_B=e(b)
matrix POOLED_V=e(V)
matrix STRATIFIED_REF=rcs_reference_V_stratified
quietly jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) preperiods(1 2) sampling(stratified)
assert "`e(rcs_sampling)'"=="stratified"
mata: assert(max(abs(st_matrix("POOLED_B")-st_matrix("e(b)")))<1e-12)
mata: assert(max(abs(st_matrix("STRATIFIED_REF")-st_matrix("e(V)")))<1e-8)
mata: assert(min(symeigenvalues(st_matrix("POOLED_V")-st_matrix("e(V)")))> -1e-12)
di "RCS_STRATIFIED_SANDWICH_REFERENCE_OK"
estimates store strat_source
jwselect, select(wbic)
assert rowsof(e(output_map))==1
estimates restore strat_source
capture noisily jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) sampling(survey)
assert _rc==198
local original_rng `c(rngstate)'
matrix BEFORE=e(b)
capture noisily jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(8) baseline(3)
assert _rc!=0
assert "`e(cmd)'"=="jwselect_rcs"
mata: assert(max(abs(st_matrix("BEFORE")-st_matrix("e(b)")))==0)
assert `"`original_rng'"'==`"`c(rngstate)'"'
capture noisily jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) preperiods(4)
assert _rc==198
capture noisily jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) pscoretrim(.995)
assert _rc==198
di "RCS_FAILURE_STATE_AND_TRIMMING_OK"
quietly jwselect_rcs y, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3)
assert e(jw_K)==0
assert colsof(e(jw_V))==1
jwselect, select(wbic)
assert rowsof(e(output_map))==1
di "RCS_ZERO_PRE_NO_COVARIATES_OK"
quietly jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(2) preperiods(1)
assert e(anticipation)==1
jwselect, select(wbic)
assert rowsof(e(output_map))==1
di "RCS_ANTICIPATION_OK"
quietly jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) preperiods(1 2)
replace x=x+1 in 1
capture noisily _jwselect_rcs
assert _rc==459
di "RCS_SIGNATURE_OK"
di "RCS_BRANCHES_ALL_OK"
log close
