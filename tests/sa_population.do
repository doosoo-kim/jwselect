version 17.0
clear all
set more off
capture log close _all
log using "tests/sa_population.log", text replace
adopath ++ "."
do _jwselect_sa_iw_mata.do
set obs 600
generate long id=_n
generate byte g=cond(_n<=300,0,cond(_n<=500,4,6))
set seed 20261104
generate double u=rnormal()
expand 8
bysort id: generate byte t=_n
generate double y=u+.1*t+(g==4 & t>=4)*1+(g==6 & t>=6)*3+rnormal()
generate byte sample=1
sort id t
mata: jw_sa_iw_inputs("y","g","sample",8,1,(0,2),J(0,1,.),0,"cohort", "D","H","V","HM","W","TW","BASES","TARGETMAP","K","SV","FV")
assert scalar(K)==6
mata:
W=st_matrix("W"); V=st_matrix("V")
assert(max(abs(W[,4]-(2/3\1/3)))<1e-12)
assert(abs(st_matrix("D")-W[,4]'*W[,5])<1e-12)
assert(abs(V[1,1]-st_numscalar("SV")-st_numscalar("FV"))<1e-12)
assert(st_numscalar("SV")>0)
/* Independent multinomial delta-method variance of estimated shares. */
p=W[,3]; w=W[,4]; d=W[,5]; psi=w'*d
vshare=sum(p:*(d:-psi):^2)/sum(p)^2/600
assert(abs(vshare-st_numscalar("SV"))<1e-12)
assert(rows(st_matrix("HM"))==6)
assert(rows(st_matrix("TW"))==6)
assert(abs(sum(st_matrix("TW")[,3])-1)<1e-12)
assert(min(symeigenvalues(V))>0)
end
display "PASS SA population: heterogeneous shares and complete joint covariance."

// A singleton supported cohort has no aggregation-share uncertainty.
mata: jw_sa_iw_inputs("y","g","sample",8,1,(0,4),J(0,1,.),0,"cohort", "D","H","V","HM","W","TW","BASES","TARGETMAP","K","SV","FV")
assert scalar(SV)==0
mata: assert(st_matrix("W")[1,1]==4 & rows(st_matrix("W"))==1)
display "PASS SA population: unsupported cohorts omitted, singleton shares fixed."

// Explicit unsupported cohort fails instead of changing the user's target.
capture mata: jw_sa_iw_inputs("y","g","sample",8,1,(0,4),6,0,"cohort", "D","H","V","HM","W","TW","BASES","TARGETMAP","K","SV","FV")
assert _rc==3498

// Noiseless equal effects give exactly zero share uncertainty.
replace y=u+.1*t+(g>0 & t>=g)*2
mata: jw_sa_iw_inputs("y","g","sample",8,1,(0,2),J(0,1,.),0,"cohort", "D","H","V","HM","W","TW","BASES","TARGETMAP","K","SV","FV")
assert scalar(SV)<1e-25
assert abs(D[1,1]-2)<1e-12
display "PASS SA population: homogeneous-effects limit."

// Full optional dependency integration: no install or update is attempted.
local absent 0
foreach cmd in eventstudyinteract reghdfe avar ftools {
    capture which `cmd'
    if _rc local absent 1
}
if `absent' {
    display "SKIP SA wrapper integration: optional dependencies absent."
}
else {
    replace y=u+.1*t+(g==4 & t>=4)*1+(g==6 & t>=6)*3+ ///
        (g==4 & t==1)*.4+rnormal()
    quietly jwselect_sa y, ivar(id) tvar(t) gvar(g)
    estimates store SA_SOURCE
    quietly jwselect_sa_iw, postwindow(0 2)
    assert "`e(cmd)'"=="jwselect_sa_iw"
    assert e(jw_nunits)==600
    assert e(iw_share_variance)>0
    assert e(N)==4800 & e(N_clusters)==600
    assert "`e(jw_contract)'"=="scalar_joint_0.9.0"
    mata: assert(rows(st_matrix("e(jw_V)"))==7)
    mata: assert(rows(st_matrix("e(jw_bases)"))==2)
    assert e(sample)==sample
    checkestimationsample
    display "PASS SA population source wrapper: signed sample and scalar contract."
    estimates store IW_SOURCE
    matrix IW_D=e(jw_d)
    matrix IW_V=e(jw_V)
    quietly jwselect, select(manual)
    assert e(scalar_target)==1
    mata: assert(cols(st_matrix("e(b)"))==1)
    mata: assert(abs(st_matrix("e(b)")[1,1]-st_matrix("IW_D")[1,1])<1e-12)
    mata: assert(abs(st_matrix("e(V)")[1,1]-st_matrix("IW_V")[1,1])<1e-12)
    quietly estat design
    mata: assert(rows(st_matrix("r(postmap)"))==1)
    capture noisily estat cells
    assert _rc==198
    assert e(scalar_target)==1
    mata: assert(abs(st_matrix("e(b)")[1,1]-st_matrix("IW_D")[1,1])<1e-12)
    display "PASS SA scalar diagnostics: design labels and cell-output refusal."
    estimates restore IW_SOURCE
    quietly jwselect, select(wbic) candidates(all)
    assert e(scalar_target)==1
    assert e(K)==6
    estimates restore IW_SOURCE
    quietly jwselect, select(focused) candidates(all) inference(none)
    assert e(scalar_target)==1
    estimates restore IW_SOURCE
    quietly jwselect, select(wbic) candidates(all) inference(wbic) wbicmethod(bootstrap) reps(1999) seed(1704)
    assert e(scalar_target)==1
    mata: assert(rows(st_matrix("e(ci)"))==1 & !hasmissing(st_matrix("e(ci)")))
    display "PASS SA population: baseline, WBIC, focused and bootstrap scalar pipeline."
    estimates restore SA_SOURCE
    quietly jwselect_sa_iw, postwindow(0 2) cohorts(4) anticipation(2)
    assert e(jw_K)==0
    assert e(iw_share_variance)==0
    quietly jwselect, select(wbic)
    assert e(K)==0 & e(scalar_target)==1
    display "PASS SA population: empty optional-restriction scalar path."
}
log close
