*! jwselect 0.5.0: actual logit/Poisson adapter checks.
*! Run from the package root in a separate Stata session. Clears memory.
*! No dependency is installed or updated by this test.
version 16.1
set more off
adopath ++ "."

local missing
foreach cmd in jwdid jwdid_estat reghdfe hdfe ftools {
    capture which `cmd'
    if _rc local missing `missing' `cmd'
}
if `"`missing'"' != "" {
    display as text "SKIP nonlinear integration; missing installed dependencies:`missing'"
    exit
}

capture program drop _jwselect_nl_metadata
program define _jwselect_nl_metadata, eclass
    syntax [, CMD(string) METHOD(string) CONVERGED(integer -1) INDEX(string)]
    if `"`cmd'"'!="" ereturn local cmd `cmd'
    if `"`method'"'!="" ereturn local cmd2 `method'
    if `converged'>=0 ereturn scalar converged=`converged'
    if `"`index'"'!="" ereturn local index `index'
end

capture mata: mata drop _jwselect_nl_canonical()
capture mata: mata drop _jwselect_nl_projection()
mata:
void _jwselect_nl_canonical(string scalar bname, string scalar vname,
                            string scalar cname)
{
    real matrix C, B, V, A
    real colvector D, H, base, ix
    real scalar g, t
    C=st_matrix(cname)
    B=st_matrix(bname)
    V=st_matrix(vname)
    D=J(0,1,.)
    H=J(0,1,.)
    for (g=3; g<=5; g=g+2) {
        base=selectindex((C[,1]:==g):&(C[,2]:==g-1))
        assert(length(base)==1)
        assert(B[1,base[1]]==0)
        assert(max(abs(V[base[1],.]))==0)
        for (t=g; t<=g+1; t++)
            D=D\selectindex((C[,1]:==g):&(C[,2]:==t))
        for (t=1; t<g-1; t++)
            H=H\selectindex((C[,1]:==g):&(C[,2]:==t))
    }
    ix=D\H
    st_matrix("JW_NL_Z",B[,ix])
    st_matrix("JW_NL_ZV",V[ix,ix])
    A=(.5,.5,0,0\0,0,.5,.5\.25,.25,.25,.25)
    st_matrix("JW_NL_A",A)
    st_matrix("JW_NL_EXPECTED_B",B[,D]*A')
    st_matrix("JW_NL_EXPECTED_V",A*V[D,D]*A')
}
void _jwselect_nl_projection()
{
    real matrix V, L, A, B, C
    V=st_matrix("JW_NL_ZV")
    A=st_matrix("JW_NL_A")
    L=(I(4),-V[1..4,5..8]*invsym(V[5..8,5..8]))
    B=A*L*st_matrix("JW_NL_Z")'
    C=A*L*V*L'*A'
    assert(max(abs(st_matrix("e(b)")'-B))<1e-8)
    assert(max(abs(st_matrix("e(V)")-C))<1e-8)
}
end

clear
set seed 681914
set obs 3600
generate long id = ceil(_n/6)
generate byte time = mod(_n-1,6)+1
generate byte group = cond(mod(id,3)==0,0,cond(mod(id,3)==1,3,5))
generate double x = sin(id/8)
generate double eta = -.8+.12*time+.15*x+.2*(group==3) ///
    +.3*(group>0 & time>=group)+.12*x*(group>0 & time>=group)
generate byte binary = runiform()<invlogit(eta)
generate double count = rpoisson(exp(eta+.5))
generate double original_order = runiform()
sort original_order

foreach link in logit poisson {
    local y binary
    if "`link'"=="poisson" local y count
    quietly jwdid `y' x, ivar(id) tvar(time) gvar(group) never ///
        method(`link') cluster(id)
    assert `"`e(cmd2)'"'=="`link'"
    assert e(converged)==1
    tempname source sourceb sourceV rawb rawV rawC condB condV
    quietly estimates store `source'
    matrix `sourceb'=e(b)
    matrix `sourceV'=e(V)
    local sourcecmd `"`e(cmdline)'"'

    capture drop __attgt__
    capture label drop __attgt__
    generate byte __attgt__=99
    label define __attgt__ 99 "nonlinear restoration sentinel"
    label values __attgt__ __attgt__
    capture drop original_esample
    generate byte original_esample=e(sample)
    quietly datasignature
    local signature `"`r(datasignature)'"'
    quietly ds
    local signature_vars `"`r(varlist)'"'
    // estimates restore can recreate/reorder Stata's own _est_* sample
    // storage. Use the same explicit variable order for the complete-exercise
    // check; the immediate adapter check below also checks variable order.
    quietly _datasignature `signature_vars'
    local ordered_signature `"`r(datasignature)'"'
    local before_rng `c(rngstate)'

    quietly _jwselect_jwdid, vce(unconditional)
    matrix `rawb'=r(bgt)
    matrix `rawV'=r(Vgt)
    matrix `rawC'=r(cells)
    assert `"`r(model)'"'=="nonlinear"
    assert `"`r(method)'"'=="`link'"
    assert `"`r(restriction_scale)'"'=="outcome_placebo"
    assert r(N_units)==600 & r(N_obs)==3600 & r(N_clusters)==600
    assert r(N_covariates)==1 & r(gap)==1 & r(anticipation)==1
    assert `"`e(cmd)'"'=="jwdid"
    assert `"`e(cmdline)'"'==`"`sourcecmd'"'
    assert e(sample)==original_esample
    mata: assert(max(abs(st_matrix("e(b)")-st_matrix("`sourceb'")))==0)
    mata: assert(max(abs(st_matrix("e(V)")-st_matrix("`sourceV'")))==0)
    quietly datasignature
    assert `"`r(datasignature)'"'==`"`signature'"'
    local label : label __attgt__ 99
    assert `"`label'"'=="nonlinear restoration sentinel"
    assert `"`c(rngstate)'"'==`"`before_rng'"'

    // Compare the complete joint covariance, including pre/post and cross-
    // cohort entries, with the upstream one-call outcome-scale margins.
    preserve
    quietly estat attgt, vce(unconditional)
    assert `"`r(agg)'"'=="attgt"
    matrix JW_NL_DIRECT_B=r(b)
    matrix JW_NL_DIRECT_V=r(V)
    restore
    quietly estimates restore `source'
    mata: assert(max(abs(st_matrix("JW_NL_DIRECT_B")-st_matrix("`rawb'")))<1e-9)
    mata: assert(max(abs(st_matrix("JW_NL_DIRECT_V")-st_matrix("`rawV'")))<1e-9)

    quietly _jwselect_jwdid, vce(conditional)
    matrix `condB'=r(bgt)
    matrix `condV'=r(Vgt)
    mata: assert(max(abs(st_matrix("`condB'")-st_matrix("`rawb'")))<1e-9)
    preserve
    quietly estat attgt
    matrix JW_NL_COND_V=r(V)
    restore
    quietly estimates restore `source'
    mata: assert(max(abs(st_matrix("JW_NL_COND_V")-st_matrix("`condV'")))<1e-9)

    // Independent point-estimate check: predict both linear indices, apply
    // the appropriate inverse link, and average the individual contrasts.
    preserve
    quietly _jwselect_nl_metadata, cmd(`link')
    quietly replace __tr__=1
    quietly predict double nl_xb1, xb
    quietly replace __tr__=0
    quietly predict double nl_xb0, xb
    if "`link'"=="logit" generate double nl_effect=invlogit(nl_xb1)-invlogit(nl_xb0)
    else generate double nl_effect=exp(nl_xb1)-exp(nl_xb0)
    forvalues j=1/`=rowsof(`rawC')' {
        local g=`rawC'[`j',1]
        local t=`rawC'[`j',2]
        quietly summarize nl_effect if group==`g' & time==`t', meanonly
        assert abs(r(mean)-`rawb'[1,`j'])<1e-8
    }
    restore
    quietly estimates restore `source'

    // Construct the canonical target/placebo vector independently. All
    // included baseline margins must already be structural zeros.
    mata: _jwselect_nl_canonical("`rawb'", "`rawV'", "`rawC'")

    quietly jwselect, select(manual) postwindow(0 1)
    assert e(K)==4 & e(P)==4 & e(G)==2
    mata: assert(max(abs(st_matrix("e(z)")-st_matrix("JW_NL_Z")))<1e-9)
    mata: assert(max(abs(st_matrix("e(V_z)")-st_matrix("JW_NL_ZV")))<1e-9)
    mata: assert(max(abs(st_matrix("e(b)")-st_matrix("JW_NL_EXPECTED_B")))<1e-9)
    mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JW_NL_EXPECTED_V")))<1e-9)
    assert e(sample)==original_esample

    // An all-moment manual adjustment agrees with the direct covariance
    // projection; this is a moment adjustment, not a nonlinear model refit.
    quietly estimates restore `source'
    matrix JW_NL_PRE=(3,1\5,1\5,2\5,3)
    quietly jwselect, select(manual) prematrix(JW_NL_PRE) postwindow(0 1)
    mata: _jwselect_nl_projection()
    quietly estat pretrend
    assert r(df)==4

    quietly estimates restore `source'
    quietly jwselect, select(wbic) wbicscope(joint) candidates(all) postwindow(0 1)
    assert `"`e(cmd)'"'=="jwselect"
    mata: assert(!hasmissing(st_matrix("e(b)")))

    // Incompatible rebasing is refused by the public adapter dispatch.
    quietly estimates restore `source'
    capture noisily jwselect, baseline(common) select(manual) postwindow(0 1)
    assert _rc==198
    assert `"`e(cmd)'"'=="jwdid" & `"`e(cmd2)'"'=="`link'"
    capture noisily jwselect, anticipation(1) select(manual) postwindow(0 1)
    assert _rc==198
    assert `"`e(cmd)'"'=="jwdid" & `"`e(cmd2)'"'=="`link'"

    // Failed validation preserves the data and source results.
    quietly _jwselect_nl_metadata, converged(0)
    capture noisily _jwselect_jwdid
    assert _rc==498
    assert e(converged)==0 & `"`e(cmd)'"'=="jwdid"
    quietly estimates restore `source'
    quietly _jwselect_nl_metadata, method(ppmlhdfe)
    capture noisily _jwselect_jwdid
    assert _rc==198
    assert `"`e(cmd2)'"'=="ppmlhdfe"
    quietly estimates restore `source'
    quietly _jwselect_nl_metadata, index(index)
    capture noisily _jwselect_jwdid
    assert _rc==198
    assert `"`e(index)'"'=="index"
    quietly estimates restore `source'
    quietly ds
    local current_vars `"`r(varlist)'"'
    local extra_vars : list current_vars - signature_vars
    local lost_vars : list signature_vars - current_vars
    assert `"`extra_vars'`lost_vars'"'==""
    quietly _datasignature `signature_vars'
    local final_signature `"`r(datasignature)'"'
    if `"`final_signature'"'!=`"`ordered_signature'"' {
        display as error "Nonlinear integration data-restoration check failed (`link')"
        display as error `"Before: `ordered_signature'"'
        display as error `"After:  `final_signature'"'
        display as error `"Original variables: `signature_vars'"'
        quietly ds
        display as error `"Current variables:  `r(varlist)'"'
    }
    assert `"`final_signature'"'==`"`ordered_signature'"'
    quietly estimates drop `source'
    display as result "jwselect nonlinear `link' integration: PASS"
}

// Positive anticipation uses the upstream baseline g-2 without rebasing.
quietly jwdid binary x, ivar(id) tvar(time) gvar(group) never ///
    method(logit) cluster(id) anticipation(2)
quietly jwselect, select(manual) anticipation(1) postwindow(0 1)
assert e(K)==2 & e(P)==4
matrix JW_NL_BASE=e(baselines)
assert JW_NL_BASE[1,2]==1 & JW_NL_BASE[2,2]==3

program drop _jwselect_nl_metadata
mata: mata drop _jwselect_nl_canonical()
mata: mata drop _jwselect_nl_projection()
display as result "jwselect nonlinear integration checks: PASS"
