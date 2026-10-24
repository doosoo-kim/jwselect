*! jwselect 0.2.0 integration checks. Run from the package root.
*! DESTRUCTIVE TO MEMORY: uses clear all; never run on unsaved user data.
version 16.1
clear all
set more off
adopath ++ "."

// Exercise first-call autoload; do not preload Mata and mask loader bugs.
capture mata: assert(findexternal("jw_stata()") != NULL)
assert _rc != 0
do tests/raw_fixture.do

// The raw fixture includes the original structural-zero reference cells.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(manual)
assert "`e(cmd)'" == "jwselect"
assert e(G) == 2 & e(P) == 7 & e(K) == 3
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_EB_BASE"))) < 1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JW_EV_BASE"))) < 1e-10)
mata: assert(max(abs(st_matrix("e(selection)"))) == 0)
assert "`e(inference)'" == "fixed_restrictions"

// Postestimation returns must copy matrices without modifying e().
estat design
mata: assert(max(abs(st_matrix("r(baselines)")-st_matrix("e(baselines)"))) == 0)
mata: assert(max(abs(st_matrix("r(postmap)")-st_matrix("e(postmap)"))) == 0)
mata: assert(max(abs(st_matrix("r(hmap)")-st_matrix("e(hmap)"))) == 0)
estat selection
mata: assert(max(abs(st_matrix("r(selection)")-st_matrix("e(selection)"))) == 0)
estat scores
mata: assert(max(abs(st_matrix("r(scores)")-st_matrix("e(scores)"))) == 0)
estat covariance
mata: assert(max(abs(st_matrix("r(V_fixed)")-st_matrix("e(V_fixed)"))) == 0)
estat cells
mata: assert(max(abs(st_matrix("r(b)")-st_matrix("e(b_cells)"))) == 0)
mata: assert(max(abs(st_matrix("r(V_fixed)")-st_matrix("e(V_cells)"))) == 0)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_EB_BASE"))) < 1e-10)
assert "`e(cmd)'" == "jwselect"

// Now compare the actually autoloaded Mata implementation to the oracle.
mata: assert(jw_version() == "0.4.0-joint1")
do tests/reference_fixture.do

// User-picked restrictions must be used in one joint MD system.
matrix JW_PRE = (3,1 \ 4,2)
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(manual) prematrix(JW_PRE)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_EB_MANUAL"))) < 1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JW_EV_MANUAL"))) < 1e-10)

// All local restrictions pass WBIC here; the oracle uses joint projection.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(wbic)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_EB_WBIC"))) < 1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JW_EV_WBIC"))) < 1e-10)
assert "`e(inference)'" == "pointwise_oracle"

// Common baseline is q-1 = 2, including the later cohort's original zero.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(manual) baseline(common)
matrix JW_CHECK = e(b)
assert abs(JW_CHECK[1,1]-.775) < 1e-10
assert abs(JW_CHECK[1,2]-(3.2/3-.06)) < 1e-10

// Anticipation shifts the baseline and excludes potentially anticipating dates.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(manual) anticipation(1)
matrix JW_CHECK = e(b)
assert abs(JW_CHECK[1,1]-(.775-.04)) < 1e-10
assert abs(JW_CHECK[1,2]-(3.2/3-.06)) < 1e-10
assert e(K) == 1

// One eligible restriction overall, owned by cohort 4 (scalar index case).
// Local candidate ownership differs from restrictions used in a joint MD fit.
foreach family in nested all {
    foreach method in wbic focused {
        jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
            select(`method') candidates(`family') anticipation(1) ///
            postwindow(0 2) focus(cohort)
        assert e(K) == 1
        matrix JW_ONEFAMILY = e(candidates)
        assert JW_ONEFAMILY[1,1] == 3 & JW_ONEFAMILY[1,3] == 0
        matrix JW_ONESEL = e(selection)
        if "`method'" == "focused" {
            // focus(cohort) fits each cohort using its own chosen restrictions.
            assert JW_ONESEL[1,1] == 0
            assert JW_ONESEL[3,1] == JW_ONESEL[2,1]
        }
        else {
            // Wald-BIC selects cohort 4's restriction in this fixture, then
            // uses the union in joint MD for every reported output.
            mata: assert(all(st_matrix("JW_ONESEL"):==1))
        }
        mata: assert(!hasmissing(st_matrix("e(b)")))
        mata: assert(!hasmissing(st_matrix("e(ci)")))
    }
}
// A single manual row belongs to cohort 4; joint MD uses it for all outputs.
matrix JW_ONEPRE = (4,1)
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(manual) anticipation(1) prematrix(JW_ONEPRE)
assert e(K) == 1
mata: assert(all(st_matrix("e(selection)"):==1))

// User-supplied cell weights preserve cohort targets and change the aggregate.
matrix JW_WEIGHTS = (3,3,1 \ 3,4,0 \ 3,5,0 \ 3,6,0 \ 4,4,0 \ 4,5,0 \ 4,6,0)
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(manual) aggweights(JW_WEIGHTS)
matrix JW_CHECK = e(b)
assert abs(JW_CHECK[1,1]-.775) < 1e-10
assert abs(JW_CHECK[1,3]-.4) < 1e-10
assert "`e(aggregate)'" == "custom"

// Zero optional preperiods exercise empty-mask handling in every module.
matrix JW_B0 = (0,.5,0,1)
matrix JW_V0 = (0,0,0,0 \ 0,.04,0,0 \ 0,0,0,0 \ 0,0,0,.09)
matrix JW_C0 = (3,2 \ 3,3 \ 4,3 \ 4,4)
jwselect, bmat(JW_B0) vmat(JW_V0) cells(JW_C0) nunits(100) ///
    select(focused) focus(both) reps(1)
assert e(K) == 0
matrix JW_CHECK = e(b)
assert abs(JW_CHECK[1,1]-.5) < 1e-10
assert abs(JW_CHECK[1,2]-1) < 1e-10
assert abs(JW_CHECK[1,3]-.75) < 1e-10
assert abs(JW_CHECK[1,4]-.75) < 1e-10

// Raw cell order is immaterial. Permute vector, covariance and metadata together.
mata:
JW_perm = (12,2,8,1,4,5,9,11,7,6,10,3)
st_matrix("JW_BP",st_matrix("JW_B")[.,JW_perm])
st_matrix("JW_VP",st_matrix("JW_V")[JW_perm,JW_perm])
st_matrix("JW_CP",st_matrix("JW_C")[JW_perm,.])
end
jwselect, bmat(JW_BP) vmat(JW_VP) cells(JW_CP) nunits(500) select(wbic)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_EB_WBIC"))) < 1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JW_EV_WBIC"))) < 1e-10)

// Focused outputs have a fixed-action covariance, not an ordinary posted e(V).
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(focused) focus(both) candidates(all) inference(none)
assert colsof(e(b)) == 4
assert "`e(inference)'" == "none"
// Check the posted VCE contract; copying an absent e(V) need not fail.
assert has_eprop("b")
assert !has_eprop("V")
mata: assert(rows(st_matrix("e(V_fixed)")) == 4)

// Default focused inference: finite intervals, bounded p-values, RNG restored.
local old_rng `c(rngstate)'
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(focused) focus(both) candidates(all) reps(19999) seed(913) optreps(999) optstarts(2) optmaxeval(80)
assert `"`c(rngstate)'"' == `"`old_rng'"'
assert "`e(inference)'" == "reselect"
matrix JW_CI = e(ci)
matrix JW_PV = e(pvalue)
mata: assert(!hasmissing(st_matrix("JW_CI")))
mata: assert(all(st_matrix("JW_CI")[.,1]:<=st_matrix("JW_CI")[.,2]))
mata: assert(all(st_matrix("JW_PV"):>=0) & all(st_matrix("JW_PV"):<=1))
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(focused) focus(both) candidates(all) reps(19999) seed(913) optreps(999) optstarts(2) optmaxeval(80)
mata: assert(max(abs(st_matrix("e(ci)")-st_matrix("JW_CI"))) == 0)

// Explicit legacy envelope remains available.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(focused) focus(cohort) candidates(all) inference(envelope)
assert "`e(inference)'" == "envelope"

// Reselection fixtures compare actual Mata draws to independent answers.
do tests/reselect_smoke.do

// Failure must preserve the previously posted estimates and RNG state.
matrix JW_SAVED_B = e(b)
matrix JW_BAD_C = JW_C
matrix JW_BAD_C[2,2] = 1
capture noisily jwselect, bmat(JW_B) vmat(JW_V) cells(JW_BAD_C) nunits(500)
assert _rc != 0
assert "`e(cmd)'" == "jwselect"
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_SAVED_B"))) == 0)
capture noisily jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(focused) candidates(all) maxjoint(1)
assert _rc != 0
capture noisily jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(manual) anticipation(-1)
assert _rc != 0
capture noisily jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(focused) reps(10)
assert _rc != 0
assert `"`c(rngstate)'"' == `"`old_rng'"'

display as result "jwselect matrix integration checks: PASS"

// Optional expanded upstream integration; never installs dependencies.
do tests/integration_jwdid.do
