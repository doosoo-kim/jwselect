version 17
clear all
set more off
capture log close _all
log using "./tests_rcs_schema.log",text replace
adopath ++ "."
set seed 819212
set obs 2400
gen x=rnormal()
gen int time=1+floor(runiform()*4)
gen int cohort=cond(runiform()<invlogit(.2*time+.2*x),4,0)
gen y=.5*x+.2*time+(cohort==4 & time==4)+rnormal()
quietly jwselect_rcs y x, tvar(time) gvar(cohort) cohort(4) target(4) baseline(3) preperiods(1 2)
estimates store source
program define rcs_mutate_schema, eclass
    args mode
    tempname A
    if "`mode'"=="signature" ereturn local datasignaturevars ""
    if inlist("`mode'","duplicate","baseline","post","unknown") {
        matrix `A'=e(jw_hmap)
        if "`mode'"=="duplicate" matrix `A'[2,2]=`A'[1,2]
        if "`mode'"=="baseline" matrix `A'[1,2]=3
        if "`mode'"=="post" matrix `A'[1,2]=4
        if "`mode'"=="unknown" matrix `A'[1,1]=9
        ereturn matrix jw_hmap=`A'
    }
    if "`mode'"=="badbaseline" {
        matrix `A'=e(jw_bases)
        matrix `A'[1,2]=-1
        ereturn matrix jw_bases=`A'
    }
    if "`mode'"=="badtarget" {
        matrix `A'=e(jw_targetmap)
        matrix `A'[1,2]=3
        ereturn matrix jw_targetmap=`A'
    }
end
foreach mode in signature duplicate baseline post unknown badbaseline badtarget {
    estimates restore source
    rcs_mutate_schema `mode'
    capture noisily jwselect, select(wbic)
    local result=_rc
    assert `result'!=0
    di "RCS_SCHEMA_GUARD_OK `mode' `result'"
}
estimates restore source
jwselect, select(wbic)
assert rowsof(e(output_map))==1
di "RCS_SCHEMA_ALL_OK"
log close
