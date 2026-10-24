*! jwselect 0.4.0: Stata 17 release integration gate.
*! Run from this package's root in a FRESH Stata session.
*! Clears data and estimates in memory. Save your work first.
version 17.0
clear all
set more off
display as text "Testing jwselect build 0.4.0-joint1"
quietly set obs 1
generate byte jw030_test_sentinel = 37
args outdir
if `"`outdir'"' == "" {
    local stamp = string(clock("`c(current_date)' `c(current_time)'", "DMYhms"), "%18.0f")
    local stamp = strtrim("`stamp'")
    local outdir "tests/release040_results_`stamp'"
}
capture mkdir `"`outdir'"'
if _rc {
    display as error "Choose a NEW test output directory with an existing parent."
    exit 602
}
adopath ++ "."
capture log close jw030tests
log using `"`outdir'/release040_tests.log"', text name(jw030tests)
display as text "Release gate: actual Stata execution, including public worker/collector APIs."
display as text "Output: `outdir'"
do tests/release030_lock_cleanup.do
do tests/release030_smoke.do `"`outdir'"'
do tests/release030_partial_cache.do
do tests/core_optimization.do
do tests/joint040_smoke.do
assert _N==1 & jw030_test_sentinel==37
display as result "jwselect 0.4.0 release tests: PASS"
display as text "This verifies implementation consistency, not global optimization or coverage."
log close jw030tests
