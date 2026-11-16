*! jwselect 0.9.1. Run after installing the package.
version 16.1
set more off
tempname B V C active comparison
matrix `B' = (.10, 0, .35, .50, .65, -.06, .04, 0, .25, .40)
matrix `C' = (3,1 \ 3,2 \ 3,3 \ 3,4 \ 3,5 \ 4,1 \ 4,2 \ 4,3 \ 4,4 \ 4,5)
matrix `active' = (1,0,1,1,1,1,1,0,1,1)'
matrix `V' = .03 * diag(`active') + .004 * `active' * `active''

// Select once within each cohort; the reported estimator is joint MD.
jwselect, bmat(`B') vmat(`V') cells(`C') nunits(500) ///
    postwindow(0 1) select(wbic) wbicscope(cohort)

// Keep the identical selected preperiods, target definitions, and weights.
// The comparison's SE columns condition on selection; they are not new
// selection-adjusted confidence intervals or tests.
estat compare
matrix `comparison'=r(table)
estat compare, cells
estat compare, output(3)

// e(b), e(V_fixed), and the selected masks remain unchanged. To request
// centered reselection inference for the fitted joint WBIC estimator, use:
// jwselect, bmat(`B') vmat(`V') cells(`C') nunits(500) ///
//     postwindow(0 1) select(wbic) wbicscope(cohort) ///
//     inference(wbic) wbicmethod(bootstrap) reps(49999) seed(12345)
