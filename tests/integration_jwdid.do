*! jwselect 0.2.0: upstream integration checks (UNEXECUTED in build environment)
*! Run from the package root. DESTRUCTIVE TO MEMORY: save user work first.
*! No dependency is installed or updated by this test.
version 16.1
set more off
adopath ++ "."

local missing
foreach cmd in jwdid jwdid_estat reghdfe hdfe ftools {
    capture which `cmd'
    if _rc local missing `missing' `cmd'
}
if `"`missing'"' != "" {
    display as text "SKIP upstream integration; missing installed dependencies:`missing'"
    exit
}

clear
set seed 491733
set obs 4800
generate long id = ceil(_n/8)
generate byte time = mod(_n-1,8)+1
generate byte group = cond(mod(id,3)==0,0,cond(mod(id,3)==1,4,6))
generate double x1 = sin(id/9)
generate double x2 = cos(id/13)
generate double y = .02*id + .15*time + .4*x1 - .2*x2 + ///
    (.7 + .3*x1)*(group>0 & time>=group) + .8*rnormal()
generate double original_order = runiform()
sort original_order

quietly jwdid y x1 x2, ivar(id) tvar(time) gvar(group) never ///
    method(regress) cluster(id)
assert `"`e(cmd2)'"' == "regress"
tempname source sourceb sourceV rawb rawV rawC conditionalb conditionalV
quietly estimates store `source'
matrix `sourceb' = e(b)
matrix `sourceV' = e(V)
local sourcecmd `"`e(cmdline)'"'
local sourceestat `"`e(estat_cmd)'"'
local sourceiv `"`e(ivar)'"'
local sourcetv `"`e(tvar)'"'
local sourcegv `"`e(gvar)'"'
local sourcegap = e(gap)
local sourceanti = e(anticipation)

// A preexisting variable and value label must survive the upstream margins code,
// which itself drops and regenerates __attgt__ and its value label.
capture drop __attgt__
capture label drop __attgt__
generate byte __attgt__ = 99
label define __attgt__ 99 "adapter restoration sentinel"
label values __attgt__ __attgt__
generate byte original_esample = e(sample)
quietly datasignature
local signature `"`r(datasignature)'"'
quietly ds
local vars `"`r(varlist)'"'

quietly _jwselect_jwdid, vce(unconditional)
matrix `rawb' = r(bgt)
matrix `rawV' = r(Vgt)
matrix `rawC' = r(cells)
assert r(N_units) == 600 & r(N_obs) == 4800
assert r(N_clusters) == 600 & r(N_covariates) == 2
assert r(gap) == 1 & r(anticipation) == 1
assert `"`r(vce)'"' == "unconditional"
assert `"`r(ivar)'"' == "id"
assert `"`r(tvar)'"' == "time"
assert `"`r(gvar)'"' == "group"
assert `"`e(cmd)'"' == "jwdid"
assert `"`e(cmdline)'"' == `"`sourcecmd'"'
assert `"`e(estat_cmd)'"' == `"`sourceestat'"'
assert `"`e(ivar)'"' == `"`sourceiv'"'
assert `"`e(tvar)'"' == `"`sourcetv'"'
assert `"`e(gvar)'"' == `"`sourcegv'"'
assert e(gap) == `sourcegap' & e(anticipation) == `sourceanti'
assert e(sample) == original_esample
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("`sourceb'"))) == 0)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("`sourceV'"))) == 0)
quietly datasignature
assert `"`r(datasignature)'"' == `"`signature'"'
quietly ds
assert `"`r(varlist)'"' == `"`vars'"'
local sentinel : label __attgt__ 99
assert `"`sentinel'"' == "adapter restoration sentinel"

// Every treated cohort and period is present, including zero normalization rows.
assert colsof(`rawb') == 16
assert rowsof(`rawC') == 16 & colsof(`rawC') == 2
assert rowsof(`rawV') == 16 & colsof(`rawV') == 16
mata:
JW_IC = st_matrix("`rawC'")
JW_IB = st_matrix("`rawb'")
JW_IV = st_matrix("`rawV'")
assert(rows(uniqrows(JW_IC)) == rows(JW_IC))
assert(all((JW_IC[,1]:==4):|(JW_IC[,1]:==6)))
for (JW_i=1; JW_i<=rows(JW_IC); JW_i++) {
    assert(JW_IC[JW_i,2]>=1 & JW_IC[JW_i,2]<=8)
    if (JW_IC[JW_i,2] == JW_IC[JW_i,1]-1) {
        assert(abs(JW_IB[1,JW_i]) < 1e-10)
        assert(max(abs(JW_IV[JW_i,.])) < 1e-10)
    }
}
JW_IA = J(3,16,0)
JW_IA[1,.] = ((JW_IC[,1]:==4):&(JW_IC[,2]:>=4))'/5
JW_IA[2,.] = ((JW_IC[,1]:==6):&(JW_IC[,2]:>=6))'/3
JW_IA[3,.] = (JW_IA[1,.]+JW_IA[2,.])/2
st_matrix("JW_UP_EXPECTED_B",(JW_IA*JW_IB')')
st_matrix("JW_UP_EXPECTED_V",JW_IA*JW_IV*JW_IA')
end

// Compare against the actual upstream one-call joint margins result.
preserve
quietly estat attgt, vce(unconditional)
assert `"`r(agg)'"' == "attgt"
matrix JW_UP_DIRECT_B = r(b)
matrix JW_UP_DIRECT_V = r(V)
restore
mata: assert(max(abs(st_matrix("JW_UP_DIRECT_B")-st_matrix("`rawb'"))) < 1e-8)
mata: assert(max(abs(st_matrix("JW_UP_DIRECT_V")-st_matrix("`rawV'"))) < 1e-8)

quietly _jwselect_jwdid, vce(conditional)
matrix `conditionalb' = r(bgt)
matrix `conditionalV' = r(Vgt)
assert `"`r(vce)'"' == "conditional"
mata: assert(max(abs(st_matrix("`conditionalb'")-st_matrix("`rawb'"))) < 1e-8)

// Manual baseline-only estimates match direct weighted cohort-time margins.
jwselect, select(manual) vce(unconditional)
assert `"`e(cmd)'"' == "jwselect"
assert e(G) == 2 & e(P) == 8 & e(K) == 6
assert e(N_units) == 600 & e(N_clusters) == 600
assert e(N_covariates) == 2 & e(N_observations) == 4800
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_UP_EXPECTED_B"))) < 1e-8)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JW_UP_EXPECTED_V"))) < 1e-8)
assert e(sample) == original_esample

quietly estimates restore `source'
jwselect, select(wbic) vce(unconditional)
assert `"`e(inference)'"' == "pointwise_oracle"
assert colsof(e(b)) == 3
mata: assert(!hasmissing(st_matrix("e(V)")))

quietly estimates restore `source'
local before_rng `c(rngstate)'
jwselect, select(focused) focus(both) vce(unconditional) ///
    candidates(nested) reps(19999) seed(7342) optreps(999) optstarts(2) optmaxeval(80)
// Initial upstream table survives posting and keeps cohort/time coordinates.
assert colsof(e(upstream_table)) == 8
assert rowsof(e(upstream_table)) == 16
mata: assert(!hasmissing(st_matrix("e(upstream_table)")[.,1..4]))
mata: assert(max(abs(st_matrix("e(upstream_table)")[.,1..2]-st_matrix("`rawC'"))) < 1e-10)
mata: assert(max(abs(st_matrix("e(upstream_table)")[.,3]-st_matrix("`rawb'")')) < 1e-10)
mata: assert(max(abs(st_matrix("e(upstream_table)")[.,4]:^2-diagonal(st_matrix("`rawV'")))) < 1e-10)
assert `"`e(inference)'"' == "reselect"
assert colsof(e(b)) == 4
mata: assert(!hasmissing(st_matrix("e(ci)")))
assert `"`c(rngstate)'"' == `"`before_rng'"'

// Default reghdfe is available only with explicitly conditional covariance.
quietly jwdid y x1 x2, ivar(id) tvar(time) gvar(group) never cluster(id)
assert `"`e(cmd2)'"' == ""
quietly _jwselect_jwdid, vce(conditional)
assert r(N_units) == 600 & r(N_clusters) == 600
assert `"`r(vce)'"' == "conditional"
capture noisily _jwselect_jwdid, vce(unconditional)
assert _rc == 198
assert `"`e(cmd)'"' == "jwdid"
assert `"`e(cmd2)'"' == ""
jwselect, select(wbic) vce(conditional)
assert `"`e(covariance)'"' == "conditional"

// Simulate a dispatcher that swallows a margins error, without replacing any
// upstream program. Only the stored estimation command's dispatch metadata is
// temporarily changed; restoring `source' restores the genuine dispatcher.
quietly estimates restore `source'
capture program drop _jwselect_test_swallow
program define _jwselect_test_swallow, rclass
    version 16.1
    syntax anything [, *]
    capture noisily error 498
    return local test_swallowed_error "yes"
end
capture program drop _jwselect_test_dispatch
program define _jwselect_test_dispatch, eclass
    version 16.1
    ereturn local estat_cmd _jwselect_test_swallow
end
quietly _jwselect_test_dispatch
quietly datasignature
local fail_signature `"`r(datasignature)'"'
capture noisily _jwselect_jwdid, vce(unconditional)
assert _rc == 498
assert `"`e(cmd)'"' == "jwdid"
assert `"`e(estat_cmd)'"' == "_jwselect_test_swallow"
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("`sourceb'"))) == 0)
quietly datasignature
assert `"`r(datasignature)'"' == `"`fail_signature'"'
quietly estimates restore `source'
assert `"`e(estat_cmd)'"' == `"`sourceestat'"'
quietly estimates drop `source'
program drop _jwselect_test_dispatch
program drop _jwselect_test_swallow

display as result "jwselect upstream integration checks: PASS"
