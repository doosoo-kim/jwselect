*! Additional tests for jwselect 0.5.0 (build 0.5.0-adapters1), revision 3, 24oct2026.
*! Extract alongside jwselect.pkg; run in a fresh Stata 17+ session.
*! Clears data and estimates. SAVE YOUR WORK FIRST.
version 17.0
args outdir scope
if "`scope'" == "" local scope all
if !inlist("`scope'", "all", "retry", "upstream") {
    display as error "Test scope must be all, retry, or upstream."
    exit 198
}
local suites 01_jwdid_joint 02_design_invariance 03_numeric_oracle 04_saved_results
if "`scope'" == "retry" local suites 01_jwdid_joint 03_numeric_oracle
if "`scope'" == "upstream" local suites 01_jwdid_joint
foreach required in jwselect.pkg jwselect.ado _jwselect_mata.do ///
    tests/joint040_fixture.do tests_extra/01_jwdid_joint.do ///
    tests_extra/02_design_invariance.do tests_extra/03_numeric_oracle.do ///
    tests_extra/numeric_oracle_fixture.do tests_extra/04_saved_results.do {
    capture confirm file "`required'"
    if _rc {
        display as error "Missing `required'. Extract this add-on into the jwselect 0.5.0 package root."
        display as error "Then cd to that directory and run: do run_extra_tests.do"
        exit 601
    }
}
clear all
set more off
set linesize 120
adopath ++ "."
if `"`outdir'"' == "" {
    local stamp = strtrim(string(clock("`c(current_date)' `c(current_time)'", "DMYhms"), "%18.0f"))
    local outdir "tests_extra/results_`stamp'"
}
capture mkdir `"`outdir'"'
if _rc {
    display as error "Choose a NEW output directory whose parent already exists."
    exit 602
}
capture log close jwextra
log using `"`outdir'/extra_tests.log"', text name(jwextra)
display as text "Additional jwselect tests: 0.5.0-adapters1; test revision 3"
display as text "Scope: `scope'; suites: `suites'"
display as text "Started: `c(current_date)' `c(current_time)'"
display as text "Working directory: `c(pwd)'"
display as text "Output directory: `outdir'"
about
which jwselect, all
which _jwselect_jwdid, all
findfile _jwselect_mata.do
display as text "Mata loader: `r(fn)'"

// Record dependency versions. No packages are installed or updated here.
local missingdeps
foreach dependency in jwdid jwdid_estat reghdfe hdfe ftools {
    capture noisily which `dependency'
    if _rc local missingdeps `missingdeps' `dependency'
}
if "`missingdeps'" != "" {
    display as error "Required commands missing: `missingdeps'"
    display as error "Tests were not run. Use the dependencies from your successful jwdid integration run."
    log close jwextra
    exit 499
}

tempname summary
file open `summary' using `"`outdir'/summary.csv"', write text
file write `summary' "suite,status,return_code" _n
local failures 0
foreach suite of local suites {
    display as text _newline "BEGIN `suite'"
    // Each suite supplies its own data and inputs; failures do not hide later results.
    capture noisily do "tests_extra/`suite'.do"
    local suite_rc = _rc
    if `suite_rc' == 0 {
        file write `summary' "`suite',PASS,0" _n
        display as result "RESULT `suite': PASS"
    }
    else {
        local failures = `failures' + 1
        file write `summary' "`suite',FAIL,`suite_rc'" _n
        display as error "RESULT `suite': FAIL (return code `suite_rc')"
    }
}
file close `summary'
display as text _newline "Finished: `c(current_date)' `c(current_time)'"
display as text "Log: `outdir'/extra_tests.log"
display as text "Summary: `outdir'/summary.csv"
if `failures' {
    display as error "jwselect 0.5.0 additional tests: FAIL (`failures' suite(s); scope `scope')"
    log close jwextra
    exit 459
}
if "`scope'" == "all" display as result "jwselect 0.5.0 additional tests: PASS"
else display as result "jwselect 0.5.0 additional retry tests: PASS (scope `scope')"
display as text "These checks assess implementation consistency; they do not establish statistical coverage."
log close jwextra
