version 17
clear all
set more off
capture log close _all
log using "./tests_rcs.log", text replace
adopath ++ "."
set seed 18042026
set obs 6000
gen double x=rnormal()
gen int time=1+floor(4*runiform())
gen double pr=invlogit(-.3+.25*x+.2*time)
gen int cohort=cond(runiform()<pr,4,0)
gen double y=1+.6*x+.3*time+.8*(cohort==4)+2*(cohort==4 & time==4)+rnormal()
jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) preperiods(1 2)
matrix list e(b)
matrix list e(V)
assert abs(_b[g4_t4]-2)<.25
assert abs(_b[g4_t1])<.25
assert abs(_b[g4_t2])<.25
assert _b[g4_t3]==0
assert e(V)[3,3]==0
do ./tests_rcs_reference.mata
mata: rcs_test_reference()
assert scalar(rcs_reference_error)<1e-8
assert scalar(rcs_reference_b_error)<1e-8
di "RCS_FULL_SANDWICH_REFERENCE_OK " scalar(rcs_reference_error)
_jwselect_rcs
assert r(N_units)==6000
matrix b=e(bgt)
matrix V=e(Vgt)
matrix C=e(cells)
estimates store source
jwselect, select(wbic)
assert rowsof(e(output_map))==1
estimates restore source
jwselect, select(focused) inference(none)
assert rowsof(e(output_map))==1
estimates restore source
matrix P=(4,1\4,2)
jwselect, select(manual) prematrix(P)
assert rowsof(e(output_map))==1
estimates restore source
jwselect, bmat(b) vmat(V) cells(C) nunits(6000) select(wbic)
di "RCS_BASIC_OK"
log close
