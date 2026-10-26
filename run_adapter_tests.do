*! jwselect 0.5.0: optional upstream integration checks.
*! Run from the extracted package root in a separate Stata 17+ session.
*! Clears data and estimates. Does not install or update any dependency.
version 17.0
clear all
set more off
adopath ++ "."
do tests/csdid_adapter.do
do tests/integration_csdid.do
do tests/integration_jwdid.do
do tests/integration_nonlinear.do
display as result "jwselect adapter test driver completed; review PASS/SKIP lines above."
