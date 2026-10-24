*! jwselect 0.4.0: joint Wald-BIC selection and joint Wald/J pretrend testing.
*! Run from the extracted package directory. Leaves the dataset intact.
version 16.1
set more off
adopath ++ "."

// A correlated two-cohort example supplied as coefficient/covariance matrices.
// Tests use these fixed values from an independent constrained-GLS calculation.
do tests/joint040_fixture.do

// The default still selects each cohort's restrictions separately.
jwselect, bmat(J40_B) vmat(J40_V) cells(J40_C) nunits(100) ///
    select(wbic) wbicscope(cohort) candidates(all)
estat selection

// Joint selection evaluates each candidate union with the FULL covariance.
// Here the separate search selects neither preperiod, while joint selects both.
jwselect, bmat(J40_B) vmat(J40_V) cells(J40_C) nunits(100) ///
    select(wbic) wbicscope(joint) candidates(all) maxjoint(4)
estat selection
estat scores

// Test every eligible preperiod jointly, relative to the chosen baselines.
// Wald and J coincide numerically under the matching covariance convention.
estat pretrend

// Test a prespecified subset: one row per cohort/date, excluding baselines.
matrix mypre = (3,1)
estat pretrend, prematrix(mypre)

// Describe the restrictions selected from the same data. The chi-square
// p-value below is nominal: it does not adjust for searching over models.
estat pretrend, selected

// Control the joint family explicitly (columns follow estat design's hmap).
// The empty candidate is inserted if absent. Candidate order breaks exact ties.
matrix myfamily = (1,0 \ 0,1 \ 1,1)
jwselect, bmat(J40_B) vmat(J40_V) cells(J40_C) nunits(100) ///
    select(wbic) wbicscope(joint) jointmatrix(myfamily) maxjoint(4)
estat pretrend

// Focused selection has a separate mask for each reported output. Specify the
// output when describing selected restrictions; all-preperiod tests need none.
jwselect, bmat(J40_B) vmat(J40_V) cells(J40_C) nunits(100) ///
    select(focused) focus(both) candidates(all) inference(none)
matrix list e(output_map)
estat pretrend, selected output(4)
