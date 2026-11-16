version 17
clear all
set more off
adopath ++ "."
log using tests/regression_090.log, text replace
do tests/raw_fixture.do
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(manual)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_EB_BASE")))<1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JW_EV_BASE")))<1e-10)
estat pretrend
estat design
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) select(wbic)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_EB_WBIC")))<1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JW_EV_WBIC")))<1e-10)
do tests/focused_ties.do
do tests/index_blocks_core.do
display "PASS 0.9 legacy covariance, selection, focused tie and index core regression checks"
log close
