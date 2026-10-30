*! jwselect 0.7.0: joint Wald-BIC selection and joint Wald/J pretrend testing.
*! Run from the extracted package directory. Leaves the dataset intact.
version 16.1
set more off
adopath ++ "."

// A correlated two-cohort example supplied as coefficient/covariance matrices.
// Fixed matrix inputs are included here; no external dataset is required.
matrix J40_B = (0.5, 0, 0.80000000000000004, 0.5, 0, 1)
matrix J40_V = (0.040000000000000001, 0, 0.0083999999999999995, 0.036000000000000004, 0, 0.022400000000000003 \ 0, 0, 0, 0, 0, 0 \ 0.0083999999999999995, 0, 0.051840000000000004, 0.0068000000000000005, 0, 0.014400000000000003 \ 0.036000000000000004, 0, 0.0068000000000000005, 0.040000000000000001, 0, 0.023200000000000002 \ 0, 0, 0, 0, 0, 0 \ 0.022400000000000003, 0, 0.014400000000000003, 0.023200000000000002, 0, 0.063760000000000011)
matrix J40_C = (3, 1 \ 3, 2 \ 3, 3 \ 4, 2 \ 4, 3 \ 4, 4)

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
