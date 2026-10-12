*! Optional actual multi-session test; run from package root in fresh sessions.
*! Save your work: this script clears memory. See release030_README.md.
version 17.0
clear all
set more off
quietly set obs 1
args stage jobdir worker workers
if `"`jobdir'"' == "" error 198
if "`workers'" == "" local workers 4
adopath ++ "."

if "`stage'" == "prepare" {
    do tests/raw_fixture.do
    local config bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
        select(focused) candidates(all) focus(both) ///
        reps(1999) optreps(199) optstarts(2) optmaxeval(80) seed(913) cachemb(64)
    jwselect, `config' inference(reselect)
    estimates store JW030_EXPECTED
    jwselect, `config' inference(none) prepare(`"`jobdir'"')
    estimates restore JW030_EXPECTED
    estimates save `"`jobdir'/test_serial_expected.ster"'
    display as result "Prepared. Open `workers' fresh Stata sessions and run their worker commands."
    exit
}
if "`stage'" == "worker" {
    if "`worker'" == "" error 198
    log using `"`jobdir'/test_worker_`worker'.log"', text replace
    jwselect_worker using `"`jobdir'"', worker(`worker') workers(`workers') cachemb(64) progress
    log close
    exit
}
if "`stage'" == "collect" {
    do tests/release030_assertions.do
    jwselect_collect using `"`jobdir'"'
    matrix JW030_ACTUAL = r(table)
    matrix JW030_ACTUAL_DIAG = r(diagnostics)
    matrix JW030_ACTUAL_EXT = r(extrema)
    matrix JW030_ACTUAL_MAP = r(output_map)
    estimates use `"`jobdir'/test_serial_expected.ster"'
    mata: jw030_equal("JW030_ACTUAL","e(table)",1e-10)
    mata: jw030_equal("JW030_ACTUAL_DIAG","e(inference_diagnostics)",1e-10)
    mata: jw030_equal("JW030_ACTUAL_EXT","e(inference_extrema)",1e-10)
    mata: jw030_equal("JW030_ACTUAL_MAP","e(output_map)",1e-10)
    display as result "Actual concurrent worker results match serial: PASS"
    exit
}
display as error "stage must be prepare, worker, or collect"
exit 198
