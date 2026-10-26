*! Additional jwselect 0.5.0 tests: actual jwdid fits, joint WBIC and Wald/MD J.
*! Test revision 3: correct _datasignature syntax for source-variable checks.
*! Stata 17.0. Save your work first: this file replaces the dataset in memory.
*! No dependencies are installed; no log is opened. Final Stata/MP 17 PASS reported by the author.
version 17.0
set more off

local XJ_missing
foreach XJ_cmd in jwselect _jwselect_jwdid jwdid jwdid_estat reghdfe hdfe ftools {
    capture which `XJ_cmd'
    if _rc local XJ_missing `XJ_missing' `XJ_cmd'
}
if `"`XJ_missing'"' != "" {
    display as error "FAIL 01_jwdid_joint: missing installed commands:`XJ_missing'"
    exit 499
}
quietly _jwselect_load
mata: assert(jw_version()=="0.5.0-adapters1")

foreach XJ_fun in close equal maps oracle selection pretrend {
    capture mata: mata drop XJ_`XJ_fun'()
}
mata:
real scalar XJ_close(real matrix a, real matrix b, real scalar tol)
{
    real colvector ix, av, bv
    if (rows(a)!=rows(b) | cols(a)!=cols(b)) return(0)
    if (!all(vec(a:>=.):==vec(b:>=.))) return(0)
    ix = vec(selectindex(vec(a:<.)))
    if (rows(ix)==0) return(1)
    av = vec(a)
    bv = vec(b)
    return(max(abs(av[ix]-bv[ix]))<=tol*max((1,max(abs(bv[ix])))))
}

void XJ_equal(string scalar aname, string scalar bname, real scalar tol)
{
    if (!XJ_close(st_matrix(aname),st_matrix(bname),tol)) {
        printf("FAIL matrix comparison: %s versus %s\n",aname,bname)
        assert(0)
    }
}

// Reconstruct the canonical design from raw cohort-time margins, independently.
void XJ_maps(string scalar basekind, real scalar delta, string scalar family)
{
    real matrix rawc, rawv, bases, hm, pm, trans, map, targets, cand
    real colvector rawb, hit, ref, gi
    real scalar g, t, bg, r, k, p, nc, a, j
    rawc = st_matrix("XJ_RC")
    rawb = st_matrix("XJ_RB")'
    rawv = st_matrix("XJ_RV")
    assert(rows(rawc)==16 & cols(rawc)==2)
    assert(rows(uniqrows(rawc))==16)
    bases = J(0,2,.)
    hm = J(0,2,.)
    pm = J(0,2,.)
    nc = 1
    for (g=4; g<=6; g=g+2) {
        bg = g-1-delta
        if (basekind=="common") bg = 3-delta
        bases = bases\(g,bg)
        a = 0
        for (t=1; t<=8; t++) {
            if (t>=g) pm = pm\(g,t)
            if (t<g-delta & t!=bg) {
                hm = hm\(g,t)
                a = a+1
            }
        }
        if (family=="all") nc = nc*2^a
        else nc = nc*(a+1)
    }
    k = rows(hm)
    p = rows(pm)
    assert(p==8 & k==6-2*delta)
    assert(XJ_close(st_matrix("e(baselines)"),bases,0))
    assert(XJ_close(st_matrix("e(hmap)"),hm,0))
    assert(XJ_close(st_matrix("e(postmap)"),pm,0))
    assert(all(hm[,2]:<hm[,1]:-delta))
    if (basekind=="common") {
        // Cohort 6 must retain eligible dates after the common baseline.
        assert(sum((hm[,1]:==6):&(hm[,2]:>3-delta))==2)
    }
    map = pm\hm
    trans = J(rows(map),rows(rawc),0)
    for (r=1; r<=rows(map); r++) {
        gi = vec(selectindex(bases[,1]:==map[r,1]))
        hit = vec(selectindex((rawc[,1]:==map[r,1]):&(rawc[,2]:==map[r,2])))
        ref = vec(selectindex((rawc[,1]:==map[r,1]):&(rawc[,2]:==bases[gi[1],2])))
        assert(rows(hit)==1 & rows(ref)==1)
        trans[r,hit[1]] = 1
        trans[r,ref[1]] = -1
    }
    assert(XJ_close(st_matrix("e(z)"),(trans*rawb)',1e-8))
    assert(XJ_close(st_matrix("e(V_z)"),trans*rawv*trans',1e-8))
    targets = J(3,p+k,0)
    targets[1,1..p] = (pm[,1]:==4)'/5
    targets[2,1..p] = (pm[,1]:==6)'/3
    targets[3,.] = (targets[1,.]+targets[2,.])/2
    assert(XJ_close(st_matrix("e(targets)"),targets,1e-12))
    cand = st_matrix("e(joint_candidates)")
    assert(rows(cand)==nc & nc<=256 & cols(cand)==k)
    assert(rows(uniqrows(cand))==nc)
    assert(all((cand:==0):|(cand:==1)))
    assert(sum(cand[1,.])==0)
    if (family=="nested") {
        for (g=4; g<=6; g=g+2) {
            gi = vec(selectindex(hm[,1]:==g))
            for (j=2; j<=rows(gi); j++) {
                assert(all(cand[,gi[j-1]]:<=cand[,gi[j]]))
            }
        }
    }
}

// Full-system covariance-weighted constrained MD residual J, computed independently of jwselect.
// Free coordinates include every post effect and every untested pre contrast.
real rowvector XJ_oracle(real colvector z, real matrix v, real rowvector mask)
{
    real scalar k, m, df, wald, mdj, pv
    real colvector restrictidx, freeidx, h, resid, beta
    real rowvector freemask
    real matrix weight, design, ident, q
    m = rows(z)
    k = cols(mask)
    df = sum(mask)
    if (df==0) return((0,0,0,.))
    restrictidx = vec(selectindex(mask')):+(m-k)
    h = z[restrictidx]
    q = v[restrictidx,restrictidx]
    assert(rank(q)==df)
    wald = h'*invsym(q)*h
    assert(rank(v)==m)
    weight = invsym(v)
    freemask = J(1,m,1)
    freemask[restrictidx] = J(1,df,0)
    freeidx = vec(selectindex(freemask'))
    ident = I(m)
    design = ident[,freeidx]
    beta = invsym(design'*weight*design)*design'*weight*z
    resid = z-design*beta
    mdj = resid'*weight*resid
    pv = chi2tail(df,wald)
    assert(XJ_close(wald,mdj,1e-7))
    return((df,wald,mdj,pv))
}

void XJ_selection()
{
    real matrix cand, v, selected, scores, jointscore, target, loading
    real colvector z, ix
    real rowvector result, chosen
    real scalar r, best, bestscore, criterion, n, k
    cand = st_matrix("e(joint_candidates)")
    v = st_matrix("e(V_z)")
    z = st_matrix("e(z)")'
    n = st_numscalar("e(N_units)")
    k = st_numscalar("e(K)")
    scores = st_matrix("e(scores)")
    jointscore = select(scores,scores[,1]:==2)
    assert(rows(jointscore)==rows(cand))
    best = 1
    bestscore = .
    for (r=1; r<=rows(cand); r++) {
        result = XJ_oracle(z,v,cand[r,.])
        criterion = result[2]-result[1]*ln(n)
        assert(XJ_close(jointscore[r,3..6],(r,result[1],result[2],criterion),1e-7))
        if (criterion<bestscore) {
            best = r
            bestscore = criterion
        }
    }
    chosen = cand[best,.]
    selected = st_matrix("e(selection)")
    assert(XJ_close(selected,J(rows(selected),1,1)*chosen,0))
    target = st_matrix("e(targets)")
    loading = target
    if (sum(chosen)>0) {
        ix = vec(selectindex(chosen')):+(rows(z)-k)
        loading[,ix] = loading[,ix]-target*v[,ix]*invsym(v[ix,ix])
    }
    assert(XJ_close(st_matrix("e(loading)"),loading,1e-7))
    assert(XJ_close(st_matrix("e(b)"),(loading*z)',1e-7))
    assert(XJ_close(st_matrix("e(V)"),loading*v*loading',1e-7))
}

void XJ_pretrend(real scalar selectedtest)
{
    real rowvector mask, wanted
    real matrix hm, savedselection
    real colvector ix
    real scalar k
    k = st_numscalar("e(K)")
    mask = J(1,k,1)
    if (selectedtest) {
        savedselection = st_matrix("e(selection)")
        mask = savedselection[1,.]
    }
    wanted = XJ_oracle(st_matrix("e(z)")',st_matrix("e(V_z)"),mask)
    assert(XJ_close(st_matrix("r(table)"),wanted,1e-7))
    assert(XJ_close(st_matrix("r(mask)"),mask,0))
    assert(st_numscalar("r(df)")==sum(mask))
    assert(st_numscalar("r(postselection)")==selectedtest)
    if (sum(mask)>0) {
        hm = st_matrix("e(hmap)")
        ix = vec(selectindex(mask'))
        assert(XJ_close(st_matrix("r(preperiods)"),hm[ix,.],0))
    }
}
end

capture program drop XJ_run_pretrend
program define XJ_run_pretrend
    version 17.0
    local XJ_rng `"`c(rngstate)'"'
    local XJ_cmdline `"`e(cmdline)'"'
    local XJ_emats b V selection scores z V_z baselines hmap postmap ///
        loading targets joint_candidates b_cells V_cells
    foreach XJ_mat of local XJ_emats {
        matrix XJ_P_`XJ_mat' = e(`XJ_mat')
    }
    quietly datasignature
    local XJ_sig `"`r(datasignature)'"'
    quietly estat pretrend
    display as text "CHECK all-test metadata: set=`r(test_set)', output=" r(output) ", df=" r(df)
    assert `"`r(test_set)'"'=="all"
    assert r(output)==0
    mata: XJ_pretrend(0)
    matrix XJ_ALL = r(table)
    quietly estat pretrend, selected
    display as text "CHECK selected-test metadata: set=`r(test_set)', output=" r(output) ", df=" r(df)
    assert `"`r(test_set)'"'=="selected"
    assert r(output)==1
    mata: XJ_pretrend(1)
    matrix XJ_SEL = r(table)
    local XJ_returnmats : r(matrices)
    if r(df)==0 {
        display as text "CHECK zero-df test has no preperiods matrix: `XJ_returnmats'"
        assert strpos(" `XJ_returnmats' "," preperiods ")==0
    }
    quietly estat pretrend, selected output(3)
    display as text "CHECK selected output(3): actual output=" r(output) ", df=" r(df)
    assert r(output)==3
    mata: XJ_pretrend(1)
    mata: XJ_equal("r(table)","XJ_SEL",1e-12)
    display as text "CHECK postestimation command/scope: cmd=`e(cmd)', scope=`e(wbic_scope)'"
    assert `"`e(cmd)'"'=="jwselect"
    assert `"`e(wbic_scope)'"'=="joint"
    display as text "CHECK postestimation command line unchanged"
    if `"`e(cmdline)'"'!=`"`XJ_cmdline'"' {
        display as text `"Expected command line: `XJ_cmdline'"'
        display as text `"Actual command line: `e(cmdline)'"'
    }
    assert `"`e(cmdline)'"'==`"`XJ_cmdline'"'
    foreach XJ_mat of local XJ_emats {
        mata: XJ_equal("e(`XJ_mat')","XJ_P_`XJ_mat'",0)
    }
    display as text "CHECK postestimation RNG unchanged: " (`"`c(rngstate)'"'==`"`XJ_rng'"')
    assert `"`c(rngstate)'"'==`"`XJ_rng'"'
    quietly datasignature
    display as text `"CHECK postestimation whole-data signature: expected=`XJ_sig'; actual=`r(datasignature)'"'
    assert `"`r(datasignature)'"'==`"`XJ_sig'"'
end

clear
set seed 617430
set obs 2400
generate long id = ceil(_n/8)
generate byte time = mod(_n-1,8)+1
generate byte group = cond(mod(id,3)==0,0,cond(mod(id,3)==1,4,6))
generate double x1 = sin(id/9)
generate double x2 = cos(id/13)
generate double y = .02*id + .15*time + .4*x1 - .2*x2 + ///
    (.7+.3*x1)*(group>0 & time>=group) + .8*rnormal()
generate double XJ_order = runiform()
sort XJ_order
tempname XJ_source
local XJ_cases = 0

foreach XJ_backend in regress reghdfe {
    local XJ_method
    if "`XJ_backend'"=="regress" local XJ_method method(regress)
    quietly jwdid y x1 x2, ivar(id) tvar(time) gvar(group) never ///
        cluster(id) `XJ_method'
    display as text "CHECK upstream backend `XJ_backend': cmd=`e(cmd)', cmd2=`e(cmd2)'"
    assert `"`e(cmd)'"'=="jwdid"
    if "`XJ_backend'"=="regress" assert `"`e(cmd2)'"'=="regress"
    else assert `"`e(cmd2)'"'==""
    quietly estimates store `XJ_source'
    matrix XJ_SOURCE_B = e(b)
    matrix XJ_SOURCE_V = e(V)
    local XJ_sourcecmd `"`e(cmdline)'"'
    local XJ_sourceestat `"`e(estat_cmd)'"'
    capture drop XJ_esample
    generate byte XJ_esample = e(sample)
    capture drop __attgt__
    capture label drop __attgt__
    generate byte __attgt__ = 99
    label define __attgt__ 99 "XJ adapter restoration sentinel"
    label values __attgt__ __attgt__
    quietly datasignature
    local XJ_sig `"`r(datasignature)'"'
    local XJ_rng `"`c(rngstate)'"'

    // Compare the same explicit source-variable list when active estimates
    // change. Stata may change its own e(sample)/stored-estimate bookkeeping.
    // Keep every current named source variable except reserved temporary
    // names and _est_* sample markers; named jwdid variables remain included.
    quietly ds
    local XJ_allvars `"`r(varlist)'"'
    local XJ_stablevars
    foreach XJ_var of local XJ_allvars {
        if substr("`XJ_var'",1,5)!="_est_" & ///
            !regexm("`XJ_var'","^__[0-9A-Za-z][0-9A-Za-z][0-9A-Za-z][0-9A-Za-z][0-9A-Za-z][0-9A-Za-z]$") {
            local XJ_stablevars `XJ_stablevars' `XJ_var'
        }
    }
    display as text "Stable source variables: `XJ_stablevars'"
    quietly _datasignature `XJ_stablevars'
    local XJ_stablesig `"`r(datasignature)'"'
    if `"`XJ_stablesig'"'=="" local XJ_stablesig `"`r(_datasignature)'"'
    display as text `"CHECK initial source-data signature is nonempty: `XJ_stablesig'"'
    assert `"`XJ_stablesig'"'!=""

    foreach XJ_vce in unconditional conditional {
        if "`XJ_backend'"=="reghdfe" & "`XJ_vce'"=="unconditional" continue
        quietly estimates restore `XJ_source'
        quietly datasignature
        local XJ_adaptersig `"`r(datasignature)'"'
        quietly _jwselect_jwdid, vce(`XJ_vce')
        matrix XJ_RB = r(bgt)
        matrix XJ_RV = r(Vgt)
        matrix XJ_RC = r(cells)
        display as text "CHECK adapter counts: N_units=" r(N_units) ", N_clusters=" r(N_clusters) ", N_obs=" r(N_obs)
        assert r(N_units)==300 & r(N_clusters)==300 & r(N_obs)==2400
        display as text "CHECK adapter design: N_covariates=" r(N_covariates) ", gap=" r(gap) ", vce=`r(vce)'"
        assert r(N_covariates)==2 & r(gap)==1
        assert `"`r(vce)'"'=="`XJ_vce'"
        display as text "CHECK adapter source command: cmd=`e(cmd)'"
        assert `"`e(cmd)'"'=="jwdid"
        display as text "CHECK adapter source command line unchanged"
        assert `"`e(cmdline)'"'==`"`XJ_sourcecmd'"'
        display as text "CHECK adapter estat dispatcher: expected=`XJ_sourceestat'; actual=`e(estat_cmd)'"
        assert `"`e(estat_cmd)'"'==`"`XJ_sourceestat'"'
        quietly count if e(sample)!=XJ_esample
        display as text "CHECK adapter sample: mismatches=" r(N)
        assert e(sample)==XJ_esample
        mata: XJ_equal("e(b)","XJ_SOURCE_B",0)
        mata: XJ_equal("e(V)","XJ_SOURCE_V",0)
        display as text "CHECK adapter RNG unchanged: " (`"`c(rngstate)'"'==`"`XJ_rng'"')
        assert `"`c(rngstate)'"'==`"`XJ_rng'"'
        quietly datasignature
        display as text `"CHECK adapter whole-data signature: expected=`XJ_adaptersig'; actual=`r(datasignature)'"'
        assert `"`r(datasignature)'"'==`"`XJ_adaptersig'"'

        foreach XJ_base in cohort common {
            foreach XJ_delta in 0 1 {
                foreach XJ_family in nested all {
                    if "`XJ_backend'"=="reghdfe" & ///
                        !("`XJ_base'"=="common" & `XJ_delta'==1 & "`XJ_family'"=="all") continue
                    local XJ_options select(wbic) wbicscope(joint) ///
                        baseline(`XJ_base') anticipation(`XJ_delta') ///
                        candidates(`XJ_family') maxcand(256) maxjoint(256)
                    display as text "Case: `XJ_backend', `XJ_vce', `XJ_base', anticipation(`XJ_delta'), `XJ_family'"
                    quietly estimates restore `XJ_source'
                    quietly jwselect, `XJ_options' vce(`XJ_vce')
                    display as text "CHECK direct fit version/scope/VCE: version=`e(version)', scope=`e(wbic_scope)', covariance=`e(covariance)'"
                    assert `"`e(version)'"'=="0.5.0"
                    assert `"`e(wbic_scope)'"'=="joint"
                    assert `"`e(covariance)'"'=="`XJ_vce'"
                    display as text "CHECK direct fit dimensions: G=" e(G) ", P=" e(P) ", K=" e(K)
                    assert e(G)==2 & e(P)==8 & e(K)==6-2*`XJ_delta'
                    display as text "CHECK direct fit units/anticipation: N_units=" e(N_units) ", anticipation=" e(anticipation)
                    assert e(N_units)==300 & e(anticipation)==`XJ_delta'
                    quietly count if e(sample)!=XJ_esample
                    display as text "CHECK direct fit sample: mismatches=" r(N)
                    assert e(sample)==XJ_esample
                    display as text "CHECK direct fit independent design/selection/MD calculations"
                    mata: XJ_maps("`XJ_base'",`XJ_delta',"`XJ_family'")
                    mata: XJ_selection()
                    local XJ_compare b V selection scores z V_z baselines hmap ///
                        postmap loading targets joint_candidates b_cells V_cells
                    foreach XJ_mat of local XJ_compare {
                        matrix XJ_D_`XJ_mat' = e(`XJ_mat')
                    }
                    display as text "CHECK direct-fit pretrend suite"
                    XJ_run_pretrend
                    matrix XJ_D_ALL = XJ_ALL
                    matrix XJ_D_SEL = XJ_SEL
                    quietly count if e(sample)!=XJ_esample
                    display as text "CHECK direct fit sample after pretrend: mismatches=" r(N)
                    assert e(sample)==XJ_esample

                    quietly jwselect, bmat(XJ_RB) vmat(XJ_RV) cells(XJ_RC) ///
                        nunits(300) `XJ_options'
                    display as text "CHECK matrix fit covariance: actual=`e(covariance)'"
                    assert `"`e(covariance)'"'=="supplied"
                    display as text "CHECK matrix/direct numerical parity"
                    foreach XJ_mat of local XJ_compare {
                        mata: XJ_equal("e(`XJ_mat')","XJ_D_`XJ_mat'",1e-8)
                    }
                    display as text "CHECK matrix-fit pretrend suite"
                    XJ_run_pretrend
                    mata: XJ_equal("XJ_ALL","XJ_D_ALL",1e-8)
                    mata: XJ_equal("XJ_SEL","XJ_D_SEL",1e-8)
                    display as text "CHECK cross-fit RNG unchanged: " (`"`c(rngstate)'"'==`"`XJ_rng'"')
                    assert `"`c(rngstate)'"'==`"`XJ_rng'"'
                    quietly datasignature
                    local XJ_wholeafter `"`r(datasignature)'"'
                    if `"`XJ_wholeafter'"'!=`"`XJ_sig'"' {
                        display as text `"NOTE whole-data signature changed across active fits: before=`XJ_sig'; after=`XJ_wholeafter'"'
                        display as text "Checking the fixed source-variable list separately from Stata sample bookkeeping."
                    }
                    quietly _datasignature `XJ_stablevars'
                    local XJ_stableafter `"`r(datasignature)'"'
                    if `"`XJ_stableafter'"'=="" local XJ_stableafter `"`r(_datasignature)'"'
                    display as text `"CHECK cross-fit source-data signature: expected=`XJ_stablesig'; actual=`XJ_stableafter'"'
                    assert `"`XJ_stableafter'"'==`"`XJ_stablesig'"'
                    local XJ_label : label __attgt__ 99
                    display as text `"CHECK sentinel value label: actual=`XJ_label'"'
                    assert `"`XJ_label'"'=="XJ adapter restoration sentinel"
                    local XJ_cases = `XJ_cases'+1
                }
            }
        }
    }
    quietly estimates drop `XJ_source'
}
display as text "CHECK completed upstream cases: actual=`XJ_cases'; expected=17"
assert `XJ_cases'==17
program drop XJ_run_pretrend
display as result "PASS 01_jwdid_joint"
