*! jwselect 0.8.8. Run after installing the package.
version 16.1
set more off

// Two cohorts, first treated in periods 3 and 4. The g-1 cells are normalized to zero.
tempname B V C active
matrix `B' = (.10, 0, .35, .50, .65, -.06, .04, 0, .25, .40)
matrix `C' = (3,1 \ 3,2 \ 3,3 \ 3,4 \ 3,5 \ 4,1 \ 4,2 \ 4,3 \ 4,4 \ 4,5)
matrix `active' = (1,0,1,1,1,1,1,0,1,1)'
matrix `V' = .03 * diag(`active') + .004 * `active' * `active''
local inputs "bmat(`B') vmat(`V') cells(`C') nunits(500) postwindow(0 1)"

// Baseline: use the maintained g-1 period alone.
jwselect, `inputs' select(manual)

// Default: nested WBIC within each cohort, then a joint MD fit.
jwselect, `inputs'
estat selection
matrix list e(ci)

// Keep the WBIC estimates and compare marginal intervals.
jwselect, `inputs' inference(wbic) wbicmethod(bootstrap) seed(12345)
estat inference
jwselect, `inputs' inference(wbic) wbicmethod(mc)
estat inference
