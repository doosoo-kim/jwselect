version 17
clear all
set more off
capture log close _all
log using "tests/build_identity.log", text replace
adopath ++ "."
do tests/raw_fixture.do
// The archive is local audit evidence, not a package dependency.
capture confirm file "../archive_preexisting_090_20261104/jwselect-0.9.0/_jwselect_mata.do"
if !_rc {
    do "../archive_preexisting_090_20261104/jwselect-0.9.0/_jwselect_mata.do"
    mata: assert(jw_version()=="0.9.0-wbicsharp1")
    capture noisily jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(manual)
    assert _rc==499
    capture noisily _jwselect_load
    assert _rc==499
    mata: mata clear
    display "ARCHIVED_BUILD_REJECTED"
}
else display "SKIP optional archived development-branch check: archive not present."
// Even an incompatible module claiming the same public version must reject.
mata:
string scalar jw_version() return("0.9.1")
end
capture noisily jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(manual)
assert _rc==499
capture noisily _jwselect_load
assert _rc==499
mata: mata clear
quietly _jwselect_load
mata: assert(jw_version()=="0.9.1" & jw_build()=="0.9.1-20261116")
quietly jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(manual)
assert "`e(version)'"=="0.9.1"
assert "`e(build)'"=="0.9.1-20261116"
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_EB_BASE")))<1e-10)
// The job directory must be new; Stata tempfile creates a unique name only.
tempfile jobdir
quietly jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(focused) focus(aggregate) inference(none) prepare("`jobdir'") ///
    reps(19999) optreps(100) optstarts(1) optmaxeval(2)
quietly _jwselect_jobguard, directory("`jobdir'") action(verify)
preserve
quietly use "`jobdir'/READY.dta", clear
assert build=="0.9.1-20261116"
quietly replace build="0.9.0"
quietly save "`jobdir'/READY.dta", replace
restore
capture noisily _jwselect_jobguard, directory("`jobdir'") action(verify)
assert _rc!=0
display "BUILD_IDENTITY_PASS: same-version code without exact build and stale batch rejected; current numerics preserved."
log close
