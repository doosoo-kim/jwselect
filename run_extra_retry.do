*! Rerun only suite 01, additional-test revision 3, 24oct2026.
*! Run from the jwselect package root in a fresh Stata 17+ session.
*! Clears data and estimates through the main driver. Save your work first.
version 17.0
args outdir
if `"`outdir'"' == "" {
    local stamp = strtrim(string(clock("`c(current_date)' `c(current_time)'", "DMYhms"), "%18.0f"))
    local outdir "tests_extra/retry_results_`stamp'"
}
do run_extra_tests.do `"`outdir'"' upstream
