*! Save/restore and postestimation checks for jwselect 0.4.0, 24oct2026.
*! Run from the package root. Replaces data with a small sentinel dataset.
version 17.0
set more off
adopath ++ "."
clear
set obs 3
generate double jwx4_sentinel = 100 + _n/7
quietly datasignature
local original_data `"`r(datasignature)'"'
local original_rng `c(rngstate)'
do tests/joint040_fixture.do

// Matrix equality that handles missing entries explicitly. A missing entry
// must occupy the same position and have the same missing-value code.
capture mata: mata drop jwx4_same()
mata:
void jwx4_same(string scalar aname, string scalar bname, real scalar tol)
{
    real matrix A, B
    real colvector ix, av, bv
    A = st_matrix(aname)
    B = st_matrix(bname)
    assert(rows(A)==rows(B) & cols(A)==cols(B))
    av = vec(A)
    bv = vec(B)
    assert(all((av:<.):==(bv:<.)))
    ix = vec(selectindex(av:>=.))
    if (rows(ix)>0) assert(all(av[ix]:==bv[ix]))
    ix = vec(selectindex(av:<.))
    if (rows(ix)>0) assert(max(abs(av[ix]-bv[ix])) <= tol*(1+max(abs(bv[ix]))))
}
end

// Exercise both WBIC scopes, fixed manual restrictions, and focused selection
// without Monte Carlo inference. The purpose is persistence, not re-testing
// the focused inference optimizer covered by the release gate.
local fit1 "select(wbic) wbicscope(joint) candidates(all)"
local fit2 "select(wbic) wbicscope(cohort) candidates(all)"
local fit3 "select(manual) prematrix(J40_PRE)"
local fit4 "select(focused) focus(both) candidates(all) inference(none)"
tempfile saved
tempname held
forvalues scenario=1/4 {
    display as text "Saved-results scenario `scenario': `fit`scenario''"
    jwselect, bmat(J40_B) vmat(J40_V) cells(J40_C) nunits(100) `fit`scenario''
    assert "`e(version)'" == "0.4.0"
    assert "`e(build)'" == "0.4.0-joint1"
    local original_matrices : e(matrices)
    foreach m of local original_matrices {
        matrix X4_`m' = e(`m')
    }
    local original_method "`e(selection_method)'"
    local original_scope "`e(wbic_scope)'"
    local original_inference "`e(inference)'"
    quietly estat pretrend
    matrix X4_ALL = r(table)
    quietly estat pretrend, selected output(1)
    matrix X4_SELECTED = r(table)
    quietly estat pretrend, prematrix(J40_PRE)
    matrix X4_CUSTOM = r(table)

    estimates store `held'
    estimates save "`saved'", replace
    // Deliberately replace the active fit before each restoration.
    quietly jwselect, bmat(J40_B) vmat(J40_V) cells(J40_C) nunits(100) select(manual)
    estimates restore `held'
    foreach m of local original_matrices {
        mata: jwx4_same("e(`m')", "X4_`m'", 0)
    }
    ereturn clear
    estimates use "`saved'"
    assert "`e(cmd)'" == "jwselect"
    assert "`e(version)'" == "0.4.0"
    assert "`e(build)'" == "0.4.0-joint1"
    assert "`e(selection_method)'" == "`original_method'"
    assert "`e(wbic_scope)'" == "`original_scope'"
    assert "`e(inference)'" == "`original_inference'"
    assert "`e(estat_cmd)'" == "jwselect_estat"
    assert "`e(predict)'" == "jwselect_p"
    assert e(N_units)==100 & e(G)==2 & e(P)==2 & e(K)==2
    local restored_matrices : e(matrices)
    local lost_matrices : list original_matrices - restored_matrices
    local extra_matrices : list restored_matrices - original_matrices
    assert "`lost_matrices'" == ""
    assert "`extra_matrices'" == ""
    foreach m of local original_matrices {
        mata: jwx4_same("e(`m')", "X4_`m'", 0)
    }

    // Replay and every documented estat subcommand must work after disk load.
    jwselect
    jwselect, nocohort celleffects
    estat design
    mata: jwx4_same("r(baselines)", "e(baselines)", 0)
    mata: jwx4_same("r(postmap)", "e(postmap)", 0)
    mata: jwx4_same("r(hmap)", "e(hmap)", 0)
    estat selection
    mata: jwx4_same("r(candidates)", "e(candidates)", 0)
    mata: jwx4_same("r(selection)", "e(selection)", 0)
    estat scores
    mata: jwx4_same("r(scores)", "e(scores)", 0)
    estat covariance
    mata: jwx4_same("r(V_fixed)", "e(V_fixed)", 0)
    estat cells
    mata: jwx4_same("r(b)", "e(b_cells)", 0)
    mata: jwx4_same("r(V_fixed)", "e(V_cells)", 0)
    matrix X4_CELLS = r(table)
    mata: st_matrix("X4_EXPECTED_CELLS", (st_matrix("e(postmap)"), st_matrix("e(b_cells)")', sqrt(diagonal(st_matrix("e(V_cells)")))))
    mata: jwx4_same("X4_CELLS", "X4_EXPECTED_CELLS", 1e-13)
    estat inference
    mata: jwx4_same("r(ci)", "e(ci)", 0)
    mata: jwx4_same("r(pvalue)", "e(pvalue)", 0)
    estat pretrend
    mata: jwx4_same("r(table)", "X4_ALL", 1e-12)
    assert r(postselection)==0 & "`r(test_set)'"=="all"
    estat pretrend, selected output(1)
    mata: jwx4_same("r(table)", "X4_SELECTED", 1e-12)
    assert r(postselection)==1 & r(output)==1
    estat pretrend, prematrix(J40_PRE)
    mata: jwx4_same("r(table)", "X4_CUSTOM", 1e-12)
    assert r(postselection)==0 & "`r(test_set)'"=="custom"

    // The explicit prediction guard must remain attached to restored results.
    capture noisily predict double jwx4_prediction
    local prediction_rc = _rc
    assert `prediction_rc' == 198
    capture confirm variable jwx4_prediction
    local variable_rc = _rc
    assert `variable_rc' == 111
    foreach m of local original_matrices {
        mata: jwx4_same("e(`m')", "X4_`m'", 0)
    }
    assert _N==3 & c(k)==1
    quietly datasignature
    assert `"`r(datasignature)'"' == `"`original_data'"'
    assert `"`c(rngstate)'"' == `"`original_rng'"'
    estimates drop `held'
}
display as result "PASS 04_saved_results"
