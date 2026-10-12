*! Run from the extracted jwselect directory. Does not clear the dataset.
version 16.1
set more off
adopath ++ "`c(pwd)'"

// A small, fully specified two-cohort vector, including zero references.
// Fixture variables are Stata matrices, so no sample dataset is needed.
do tests/raw_fixture.do

// Default baseline g-1, manually retain the baseline only.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(manual)
estat design

// Select eligible preperiods separately, then estimate using their joint union.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(wbic)
estat selection

// Cohort-focused effects AND separately selected aggregate-focused effect.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(focused) focus(both) candidates(all) seed(12345)
estat covariance
estat inference

// Replay can hide cohort rows without deleting those stored results.
jwselect, nocohort
matrix list e(ci)
matrix list e(pvalue)

// Use the common baseline q-1 and a shared event-time window.
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) ///
    select(wbic) baseline(common) postwindow(0 2)
