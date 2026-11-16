*! Release contract: full conditional-index blocks are the nonlinear default.
*! Run from the package root; uses installed jwdid dependencies only.
version 17.0
clear all
set more off
capture log close _all
log using "./tests/nonlinear_index_release.log", text replace
adopath ++ "."
set seed 184903
set obs 3600
generate long id=ceil(_n/6)
generate byte time=mod(_n-1,6)+1
generate byte group=cond(mod(id,3)==0,0,cond(mod(id,3)==1,3,5))
generate double x=sin(id/8)+.2*(group==3)
generate double z=cos(id/11)-.1*(group==5)
generate double eta=-.5+.08*time+.2*x-.1*z+.2*(group==3)+ ///
    (.4+.15*x)*(group>0 & time>=group)
generate byte y=runiform()<invlogit(eta)
generate double yp=rpoisson(exp(eta))
foreach link in logit poisson {
    local yy y
    if "`link'"=="poisson" local yy yp
    foreach basis in centered raw {
        local xo
        if "`basis'"=="raw" local xo xasis
        quietly jwdid `yy' x z, ivar(id) tvar(time) gvar(group) ///
            never method(`link') cluster(id) `xo'
        estimates store NL_SOURCE
        matrix SOURCE_B=e(b)
        matrix SOURCE_V=e(V)
        local sourcecmd `"`e(cmdline)'"'
        local sourcerng `"`c(rngstate)'"'
        quietly datasignature
        local sourcesig `"`r(datasignature)'"'

        quietly _jwselect_index, vce(conditional)
        assert r(N_index)==12
        assert r(covariance_error)<2e-5
        assert r(margins_error)<5e-5
        assert "`r(restriction_scale)'"=="conditional_index"
        matrix COND_B=r(bgt)
        matrix COND_V=r(Vgt)
        matrix COND_H=r(h)
        matrix COND_VH=r(Vh)
        quietly _jwselect_index, vce(unconditional)
        mata: assert(mreldif(st_matrix("r(bgt)"),st_matrix("COND_B"))<1e-12)
        mata: assert(mreldif(st_matrix("r(h)"),st_matrix("COND_H"))<1e-12)
        mata: assert(mreldif(st_matrix("r(Vh)"),st_matrix("COND_VH"))<1e-12)
        mata: assert(mreldif(st_matrix("r(Vgt)"),st_matrix("COND_V"))>1e-8)
        mata: assert(mreldif(st_matrix("e(b)"),st_matrix("SOURCE_B"))==0)
        mata: assert(mreldif(st_matrix("e(V)"),st_matrix("SOURCE_V"))==0)
        assert `"`e(cmdline)'"'==`"`sourcecmd'"'
        assert `"`c(rngstate)'"'==`"`sourcerng'"'
        quietly datasignature
        assert `"`r(datasignature)'"'==`"`sourcesig'"'

        quietly jwselect, select(wbic) candidates(all) postwindow(0 1) inference(none)
        assert "`e(restrictions)'"=="index"
        assert "`e(restriction_scale)'"=="conditional_index"
        assert e(K)==12 & e(K_blocks)==4
        matrix DEFAULT_B=e(b)
        matrix DEFAULT_V=e(V)
        matrix DEFAULT_Z=e(z)
        matrix DEFAULT_ZV=e(V_z)
        matrix DEFAULT_S=e(selection)
        quietly estat pretrend
        assert r(df)==12
        quietly estimates restore NL_SOURCE
        quietly jwselect, restrictions(index) select(wbic) candidates(all) ///
            postwindow(0 1) inference(none)
        mata: assert(mreldif(st_matrix("e(b)"),st_matrix("DEFAULT_B"))==0)
        mata: assert(mreldif(st_matrix("e(V)"),st_matrix("DEFAULT_V"))==0)
        mata: assert(mreldif(st_matrix("e(z)"),st_matrix("DEFAULT_Z"))==0)
        mata: assert(mreldif(st_matrix("e(V_z)"),st_matrix("DEFAULT_ZV"))==0)
        mata: assert(mreldif(st_matrix("e(selection)"),st_matrix("DEFAULT_S"))==0)

        quietly estimates restore NL_SOURCE
        quietly jwselect, select(focused) focus(aggregate) candidates(all) ///
            postwindow(0 1) inference(none)
        assert e(K)==12 & e(K_blocks)==4
        matrix FOCUS_B=e(b)
        matrix FOCUS_V=e(V_fixed)
        matrix FOCUS_S=e(selection)
        quietly estimates restore NL_SOURCE
        quietly jwselect, restrictions(index) select(focused) focus(aggregate) ///
            candidates(all) postwindow(0 1) inference(none)
        mata: assert(mreldif(st_matrix("e(b)"),st_matrix("FOCUS_B"))==0)
        mata: assert(mreldif(st_matrix("e(V_fixed)"),st_matrix("FOCUS_V"))==0)
        mata: assert(mreldif(st_matrix("e(selection)"),st_matrix("FOCUS_S"))==0)
        if "`basis'"=="centered" {
            matrix CENTER_B=e(b)
            matrix CENTER_V=e(V_fixed)
            matrix CENTER_H=e(index_b)
            matrix CENTER_VH=e(V_index)
            matrix CENTER_C=e(Cov_att_index)
        }
        else {
            mata: assert(mreldif(st_matrix("e(b)"),st_matrix("CENTER_B"))<2e-5)
            mata: assert(mreldif(st_matrix("e(V_fixed)"),st_matrix("CENTER_V"))<2e-5)
            mata: assert(mreldif(st_matrix("e(index_b)"),st_matrix("CENTER_H"))<2e-5)
            mata: assert(mreldif(st_matrix("e(V_index)"),st_matrix("CENTER_VH"))<2e-5)
            mata: assert(mreldif(st_matrix("e(Cov_att_index)"),st_matrix("CENTER_C"))<2e-5)
        }
        display "PASS full-block source/default/focused covariance: `link' `basis'."
    }
}

// Cancellation distinguishes the supported full-block null from the
// explicitly requested legacy placebo diagnostic; it is not a new theory.
clear
set seed 719451
set obs 600
generate long id=_n
generate byte group=cond(id<=300,3,0)
generate double x=mod(id,3)-1
expand 4
bysort id: generate byte time=_n
generate double probability=.5
replace probability=invlogit(x*ln(3)) if group==3 & time==1
replace probability=.75 if group==3 & time>=3
generate double orderdraw=runiform()
bysort group x time (orderdraw): generate byte y=_n<=round(_N*probability)
sort id time
quietly jwdid y x, ivar(id) tvar(time) gvar(group) never method(logit) cluster(id)
estimates store CANCEL_SOURCE
quietly jwselect, select(wbic) candidates(all) inference(none)
assert e(K)==2 & e(K_blocks)==1
assert "`e(restrictions)'"=="index"
quietly estat pretrend
assert r(df)==2 & r(p)<.01
quietly estimates restore CANCEL_SOURCE
quietly jwselect, restrictions(placebo) select(wbic) candidates(all) inference(none)
assert e(K)==1 & e(K_blocks)==1
assert "`e(restrictions)'"=="placebo"
assert "`e(restriction_scale)'"=="outcome_placebo"
quietly estat pretrend
assert r(df)==1 & r(p)>.99
display "PASS full-block default: cancelling conditional departures remain testable."
display "NONLINEAR_INDEX_RELEASE_PASS"
log close
