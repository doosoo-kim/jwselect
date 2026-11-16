*! jwselect 0.9.0: real csdid 2.x trimming guard and unchanged-data tests
version 17.0
clear all
set more off
capture log close _all
log using "tests/dr_untrimmed.log", text replace
adopath ++ "."
capture which csdid
if _rc {
    display "SKIP DR untrimmed: optional csdid is not installed."
    log close
    exit
}
capture quietly csdid version
if _rc {
    display "SKIP DR untrimmed: requires optional csdid 2.x."
    log close
    exit
}
set seed 261104
set obs 4800
generate long id=ceil(_n/8)
generate byte time=mod(_n-1,8)+1
generate byte group=cond(mod(id,3)==0,0,cond(mod(id,3)==1,4,6))
generate double x=sin(id/9)
generate double y=.02*id+.15*time+.4*x+(.7+.3*x)*(group>0 & time>=group)+.8*rnormal()

// Even a seemingly benign .995 fit is outside the strict untrimmed contract.
quietly csdid y x, ivar(id) time(time) gvar(group) method(dr) ///
    nevertreated base_period(universal) analytical pointwise
assert abs(e(pscoretrim)-.995)<1e-12
matrix SAVED_B=e(b)
matrix SAVED_V=e(V)
generate byte original_sample=e(sample)
quietly datasignature
local sig `"`r(datasignature)'"'
local oldrng `c(rngstate)'
local oldcmd `"`e(cmdline)'"'
capture noisily _jwselect_csdid
assert _rc==198
assert "`e(cmd)'"=="csdid"
assert `"`e(cmdline)'"'==`"`oldcmd'"'
assert abs(e(pscoretrim)-.995)<1e-12
assert e(sample)==original_sample
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("SAVED_B")))==0)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("SAVED_V")))==0)
quietly datasignature
assert `"`r(datasignature)'"'==`"`sig'"'
assert `"`c(rngstate)'"'==`"`oldrng'"'
capture noisily jwselect, select(wbic) postwindow(0 2)
assert _rc==198
assert "`e(cmd)'"=="csdid"
assert e(sample)==original_sample
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("SAVED_B")))==0)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("SAVED_V")))==0)
quietly datasignature
assert `"`r(datasignature)'"'==`"`sig'"'
assert `"`c(rngstate)'"'==`"`oldrng'"'
display "PASS DR trim guard: default .995 rejected; source results, sample, data and RNG preserved."

// Actual untrimmed nuisance fits and full influence covariance are accepted.
quietly csdid y x, ivar(id) time(time) gvar(group) method(dr) ///
    nevertreated base_period(universal) analytical pointwise pscoretrim(1)
assert e(pscoretrim)==1
estimates store UNTRIMMED_SOURCE
quietly _jwselect_csdid
assert r(pscoretrim)==1
matrix B=r(bgt)
matrix V=r(Vgt)
matrix C=r(cells)
assert r(N_units)==600 & r(N_obs)==4800
mata:
cells=st_matrix("C")
weights=J(3,rows(cells),0)
weights[1,.]=((cells[,1]:==4):&(cells[,2]:>=4):&(cells[,2]:<=6))'/3
weights[2,.]=((cells[,1]:==6):&(cells[,2]:>=6):&(cells[,2]:<=8))'/3
weights[3,.]=(weights[1,.]+weights[2,.])/2
st_matrix("EXPECTED_B",(weights*st_matrix("B")')')
st_matrix("EXPECTED_V",weights*st_matrix("V")*weights')
end
quietly jwselect, select(manual) postwindow(0 2)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("EXPECTED_B")))<1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("EXPECTED_V")))<1e-10)
assert e(N_units)==600 & e(N_observations)==4800
assert e(sample)==original_sample
estimates restore UNTRIMMED_SOURCE
quietly jwselect, select(wbic) wbicscope(joint) candidates(all) postwindow(0 2)
matrix SELECTED_B=e(b)
matrix SELECTED_V=e(V)
quietly jwselect, bmat(B) vmat(V) cells(C) nunits(600) ///
    select(wbic) wbicscope(joint) candidates(all) postwindow(0 2)
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("SELECTED_B")))<1e-10)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("SELECTED_V")))<1e-10)
estimates restore UNTRIMMED_SOURCE
quietly jwselect, select(focused) focus(aggregate) candidates(all) postwindow(0 2) inference(none)
assert e(K)>0
display "PASS DR untrimmed: real source accepted, baseline covariance and joint WBIC matrix equivalence verified."
log close
