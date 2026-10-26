*! jwselect 0.5.0: integration with a real, separately installed csdid 2.x
*! Run from package root in a fresh test session. Installs nothing.
version 16.1
set more off
adopath ++ "."
capture which csdid
if _rc {
    display as text "SKIP CS-DR integration: csdid is not installed"
    exit
}
clear
// Legacy SSC csdid has no bare version subcommand. Check before using 2.x
// options so an older optional dependency produces an explicit skip.
capture quietly csdid version
if _rc {
    display as text "SKIP CS-DR integration: installed csdid does not expose the 2.x version interface"
    exit
}
set seed 260526
set obs 4800
generate long id = ceil(_n/8)
generate byte time = mod(_n-1,8)+1
generate byte group = cond(mod(id,3)==0,0,cond(mod(id,3)==1,4,6))
generate double x = sin(id/9)
generate double y = .02*id + .15*time + .4*x + ///
    (.7+.3*x)*(group>0 & time>=group) + .8*rnormal()
quietly csdid y x, ivar(id) time(time) gvar(group) method(dr) ///
    nevertreated base_period(universal) analytical pointwise
if substr("`e(version)'",1,2)!="2." {
    display as text "SKIP CS-DR integration: requires csdid 2.0.1 or compatible 2.x"
    exit
}
tempname SRC B V C JB JV JB2 JV2 A
quietly estimates store `SRC'
generate byte original_sample = e(sample)
quietly datasignature
local sig `"`r(datasignature)'"'
quietly _jwselect_csdid
matrix `B'=r(bgt)
matrix `V'=r(Vgt)
matrix `C'=r(cells)
assert r(N_units)==600 & r(N_obs)==4800
assert r(N_covariates)==1 & r(anticipation)==0
quietly datasignature
assert `"`r(datasignature)'"'==`"`sig'"'
assert "`e(cmd)'"=="csdid"

// A fixed, equal-cohort target is evaluated directly from unaggregated ATTs.
// It need not equal csdid's default aggregation with estimated cohort shares.
mata:
csC=st_matrix("`C'")
csW=J(3,rows(csC),0)
csW[1,.]=((csC[,1]:==4):&(csC[,2]:>=4):&(csC[,2]:<=6))'/3
csW[2,.]=((csC[,1]:==6):&(csC[,2]:>=6):&(csC[,2]:<=8))'/3
csW[3,.]=(csW[1,.]+csW[2,.])/2
st_matrix("CS_EXPECT_B",(csW*st_matrix("`B'")')')
st_matrix("CS_EXPECT_V",csW*st_matrix("`V'")*csW')
end
quietly jwselect, select(manual) postwindow(0 2)
assert "`e(source_estimator)'"=="csdid"
assert e(N_units)==600 & e(N_observations)==4800
assert e(sample)==original_sample
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("CS_EXPECT_B")))<1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("CS_EXPECT_V")))<1e-10)

// Adapter and explicit matrix input must agree under joint WBIC and manual
// restrictions, using all cross-cohort covariance terms in both paths.
quietly estimates restore `SRC'
quietly jwselect, select(wbic) wbicscope(joint) candidates(all) postwindow(0 2)
matrix `JB'=e(b)
matrix `JV'=e(V)
quietly estat pretrend
scalar CS_W=r(chi2)
quietly jwselect, bmat(`B') vmat(`V') cells(`C') nunits(600) ///
    select(wbic) wbicscope(joint) candidates(all) postwindow(0 2)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("`JB'")))<1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("`JV'")))<1e-10)
quietly estat pretrend
assert abs(r(chi2)-CS_W)<1e-10

// A direct csdid pretrend Wald uses the same preperiod estimates/covariance.
quietly estimates restore `SRC'
scalar CS_UP_W=e(wald_stat)
quietly jwselect, select(manual) postwindow(0 2)
quietly estat pretrend
assert abs(r(chi2)-CS_UP_W)<1e-8

quietly estimates restore `SRC'
capture noisily jwselect, baseline(common)
assert _rc==198
assert "`e(cmd)'"=="csdid"
capture noisily jwselect, anticipation(1)
assert _rc==198
assert "`e(cmd)'"=="csdid"

// Active aggregation posting is rejected. Plain estat without post is safe.
quietly estat event
quietly _jwselect_csdid
quietly estimates restore `SRC'
quietly estat event, post
capture noisily _jwselect_csdid
assert _rc==498
quietly estimates restore `SRC'

// Explicit ivar clustering retains the whole covariance and source sample.
quietly csdid y x, ivar(id) time(time) gvar(group) method(dr) ///
    nevertreated base_period(universal) analytical pointwise cluster(id)
quietly _jwselect_csdid
assert r(N_clusters)==600

// An anticipated-treatment fit has a different, verified normalization.
quietly csdid y x, ivar(id) time(time) gvar(group) method(dr) ///
    nevertreated base_period(universal) analytical pointwise anticipation(1)
quietly _jwselect_csdid
assert r(anticipation)==1
matrix `B'=r(bgt)
matrix `V'=r(Vgt)
matrix `C'=r(cells)
quietly jwselect, select(manual) anticipation(1) postwindow(0 2)
assert e(anticipation)==1

// Changed defaults are never silently accepted as the intended design.
quietly csdid y x, ivar(id) time(time) gvar(group) method(dr) ///
    base_period(universal) analytical pointwise
capture noisily _jwselect_csdid
assert _rc==198
quietly csdid y x, ivar(id) time(time) gvar(group) method(dr) ///
    nevertreated base_period(varying) analytical pointwise
capture noisily _jwselect_csdid
assert _rc==198
quietly csdid y x, ivar(id) time(time) gvar(group) method(dr) ///
    nevertreated base_period(universal) pointwise reps(25) rseed(260526)
capture noisily _jwselect_csdid
assert _rc==198

quietly estimates drop `SRC'
display as result "jwselect 0.5.0 actual CS-DR integration checks: PASS"
