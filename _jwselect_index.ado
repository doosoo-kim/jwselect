*! jwselect 0.7.0 30oct2026: complete nonlinear conditional-index blocks
program define _jwselect_index, rclass
    version 15.1
    syntax [, VCE(string)]
    if `"`vce'"' == "" local vce unconditional
    local vce = lower(`"`vce'"')
    if !inlist(`"`vce'"', "conditional", "unconditional") {
        di as err "vce() must be conditional or unconditional"
        exit 198
    }
    if `"`e(cmd)'"' != "jwdid" {
        di as err "The index adapter requires active jwdid results"
        exit 301
    }
    if !inlist(`"`e(cmd2)'"', "logit", "poisson") {
        di as err "restrictions(index) requires jwdid with method(logit) or method(poisson)"
        exit 198
    }

    // First reuse the outcome adapter's complete specification, support, and
    // sample checks. Its conditional margins provide an independent numerical
    // Jacobian check. The joint influence covariance below is authoritative.
    tempname held B V C T H VH VDH HM IT
    quietly estimates store `held'
    capture noisily _jwselect_jwdid, vce(conditional)
    local rc = _rc
    if !`rc' {
        matrix `B' = r(bgt)
        matrix `V' = r(Vgt)
        matrix `C' = r(cells)
        matrix `T' = r(upstream_table)
        foreach key in N_units N_obs N_clusters N_covariates gap anticipation df_r {
            local `key' = r(`key')
        }
        foreach key in ivar tvar gvar sourcecmd sourcevce upstream_stat model method {
            local `key' `"`r(`key')'"'
        }
        preserve
        capture noisily _jwselect_index_work, vce(`vce') bname(`B') ///
            vname(`V') cname(`C') hname(`H') vhname(`VH') ///
            vdhname(`VDH') hmname(`HM') itname(`IT')
        local rc = _rc
        if !`rc' {
            foreach key in N_index covariance_factor covariance_error margins_error {
                local `key' = r(`key')
            }
            local index_covariates `"`r(index_covariates)'"'
        }
        restore
    }
    quietly estimates restore `held'
    quietly estimates drop `held'
    if `rc' exit `rc'

    return clear
    return matrix bgt = `B'
    return matrix Vgt = `V'
    return matrix cells = `C'
    // This is the actual upstream conditional comparison table, not a table
    // relabelled with our unconditional influence-function standard errors.
    return matrix upstream_table = `T'
    if `N_index' > 0 {
        return matrix h = `H'
        return matrix Vh = `VH'
        return matrix Vdh = `VDH'
        return matrix hmap = `HM'
        return matrix index_terms = `IT'
    }
    foreach key in N_units N_obs N_clusters N_covariates gap anticipation df_r ///
        N_index covariance_factor covariance_error margins_error {
        return scalar `key' = ``key''
    }
    foreach key in ivar tvar gvar sourcecmd sourcevce upstream_stat model method {
        return local `key' `"``key''"'
    }
    return local index_covariates `"`index_covariates'"'
    return local restriction_scale conditional_index
    return local index_basis raw_covariates
    return local upstream_vce conditional
    return local covariance_method stacked_unit_influence
    return local vce `"`vce'"'
end

program define _jwselect_index_work, rclass
    version 15.1
    syntax , VCE(string) BNAME(name) VNAME(name) CNAME(name) ///
        HNAME(name) VHNAME(name) VDHNAME(name) HMNAME(name) ITNAME(name)
    local backend `"`e(cmd2)'"'
    local ivar `"`e(ivar)'"'
    local tvar `"`e(tvar)'"'
    local gvar `"`e(gvar)'"'
    local gap = e(gap)
    local anti = e(anticipation)
    local sourcecmd `"`e(cmdline)'"'
    gettoken cmd commandargs : sourcecmd
    quietly _jwselect_index_parse `commandargs'
    local xs `"`r(covariates)'"'
    local y `"`r(outcome)'"'

    tempname BB VV OO MB nindex factor coverror margerror
    matrix `BB' = e(b)
    matrix `VV' = e(V)
    capture matrix `MB' = e(V_modelbased)
    if _rc {
        di as err "The nonlinear fit does not expose the model-based Hessian covariance"
        exit 498
    }
    quietly _ms_omit_info `BB'
    matrix `OO' = r(omit)
    tempvar sample uid xb one
    quietly generate byte `sample' = e(sample)
    quietly keep if `sample'
    quietly egen long `uid' = group(`ivar')
    quietly sort `uid' `tvar'
    quietly generate double `one' = 1
    quietly _predict double `xb', xb
    // In never mode the current upstream release sets __tr__=1 for all
    // observations; only treated-cohort interactions enter the design.
    capture assert __tr__ == 1 if `gvar'>0
    if _rc {
        di as err "jwdid's treatment-index variable has changed or uses an unsupported normalization"
        exit 498
    }

    local names : colnames `BB'
    local design
    local active
    forvalues j=1/`=colsof(`BB')' {
        if `OO'[1,`j'] == 0 {
            local name : word `j' of `names'
            local active `active' `j'
            if `"`name'"' == "_cons" local design `design' `one'
            else {
                quietly fvrevar `name'
                local fv `"`r(varlist)'"'
                if `: word count `fv'' != 1 {
                    di as err "Cannot reconstruct identified model column: `name'"
                    exit 498
                }
                tempvar dx
                quietly generate double `dx' = `fv'
                local design `design' `dx'
            }
        }
    }
    // All Mata state is local to this call. No data, fit, or RNG is modified
    // outside the outer preserve/estimate hold, including on validation errors.
    mata: _jwidx_build("`BB'", "`VV'", "`MB'", "`OO'", ///
        "`bname'", "`vname'", "`cname'", "`hname'", "`vhname'", ///
        "`vdhname'", "`hmname'", "`itname'", "`nindex'", "`factor'", ///
        "`coverror'", "`margerror'", "`design'", "`active'", ///
        "`uid'", "`y'", "`xb'", "`gvar'", "`tvar'", "`xs'", ///
        "`backend'", "`vce'", `gap', `anti')
    if scalar(`nindex') > 0 {
        local hnames
        forvalues j=1/`=scalar(`nindex')' {
            local hnames `hnames' h`j'
        }
        matrix colnames `hname' = `hnames'
        matrix colnames `vhname' = `hnames'
        matrix rownames `vhname' = `hnames'
        matrix colnames `vdhname' = `hnames'
        matrix rownames `hmname' = `hnames'
        matrix colnames `hmname' = cohort time
        matrix rownames `itname' = `hnames'
        matrix colnames `itname' = source_coef component
    }
    return scalar N_index = scalar(`nindex')
    return scalar covariance_factor = scalar(`factor')
    return scalar covariance_error = scalar(`coverror')
    return scalar margins_error = scalar(`margerror')
    return local index_covariates `"`xs'"'
end

// The outcome adapter already validated every option. Its private subprograms
// are not callable from another ado file, so extract only the plain variables.
program define _jwselect_index_parse, rclass
    version 15.1
    syntax varlist(numeric) [if] [in] [, *]
    gettoken y xs : varlist
    return local outcome `"`y'"'
    return local covariates `"`xs'"'
end

mata:
void _jwidx_fail(string scalar message)
{
    errprintf("%s\n", message)
    _error(498)
}

real colvector _jwidx_mean(real colvector eta, string scalar method)
{
    if (method == "logit") return(invlogit(eta))
    return(exp(eta))
}

real colvector _jwidx_deriv(real colvector mu, string scalar method)
{
    if (method == "logit") return(mu :* (1 :- mu))
    return(mu)
}

// Each active treatment coefficient maps to (cohort, time, component,
// centered). Structural b/o nuisance columns have already been removed.
real matrix _jwidx_terms(string colvector names, string scalar gvar,
    string scalar tvar, string rowvector xs)
{
    real matrix ans
    real scalar j, k, n, g, t, comp, cent, found, tr
    string rowvector parts
    string scalar part
    ans = J(rows(names), 4, .)
    for (j=1; j<=rows(names); j++) {
        if (!strpos(names[j], "__tr__")) continue
        parts = tokens(subinstr(names[j], "#", " "))
        g = t = .
        comp = 0
        cent = 0
        tr = 0
        for (k=1; k<=cols(parts); k++) {
            part = parts[k]
            if (part == "c.__tr__") {
                tr++
                continue
            }
            if (regexm(part, "^([0-9]+)[bon]*[.]"+gvar+"$")) {
                g = strtoreal(regexs(1))
                continue
            }
            if (regexm(part, "^([0-9]+)[bon]*[.]"+tvar+"$")) {
                t = strtoreal(regexs(1))
                continue
            }
            found = 0
            for (n=1; n<=cols(xs); n++) {
                if (part == "c."+xs[n] | part == "c._x_"+xs[n]) {
                    if (comp != 0) _jwidx_fail("Unexpected higher-order treatment interaction")
                    comp = n
                    cent = (part == "c._x_"+xs[n])
                    found++
                }
            }
            if (found != 1) _jwidx_fail("Unrecognized treatment-index coefficient: "+names[j])
        }
        if (tr != 1 | missing(g) | missing(t) | g<=0)
            _jwidx_fail("Incomplete treatment-index coefficient: "+names[j])
        ans[j,] = (g,t,comp,cent)
    }
    return(ans)
}

void _jwidx_build(string scalar bbname, string scalar vvname,
    string scalar mbname, string scalar ooname, string scalar bname,
    string scalar vname, string scalar cname, string scalar hname,
    string scalar vhname, string scalar vdhname, string scalar hmname,
    string scalar itname, string scalar nkname, string scalar facname,
    string scalar errname, string scalar merrname, string scalar design,
    string scalar actives, string scalar uid, string scalar yvar,
    string scalar xbvar, string scalar gvar, string scalar tvar,
    string scalar xvars, string scalar method, string scalar vce,
    real scalar gap, real scalar anti)
{
    real matrix X, X0, XS, Z, H, Hinv, Q, V0, Vs, MB, panel, term,
        cells, A, D, L, HM, IT, QH, QD, joint, VV, Vmarg, tmp,
        bfull, vfull, om
    real colvector y, id, g, t, eta, mu, w, eta0, mu0, w0, idx,
        selectrow, coeff, f, contrib, lead, gc, src, scoresum
    real rowvector active, beta, b, meanx, pick, h
    real scalar n, p, m, q, N, j, k, c, r, gval, tval, ng,
        nlead, K, fac, err, margerr, denom, val, iscenter, zcol
    string matrix stripe
    string rowvector xs

    active = strtoreal(tokens(actives))
    bfull = st_matrix(bbname)
    vfull = st_matrix(vvname)
    om = st_matrix(ooname)
    beta = bfull[,active]
    Vs = vfull[active,active]
    MB = st_matrix(mbname)[active,active]
    X = st_data(.,tokens(design))
    y = st_data(.,yvar)
    id = st_data(.,uid)
    g = st_data(.,gvar)
    t = st_data(.,tvar)
    xs = tokens(xvars)
    q = cols(xs)
    XS = J(rows(X),0,.)
    if (q) XS = st_data(.,xs)
    n = rows(X)
    p = cols(X)
    cells = st_matrix(cname)
    b = st_matrix(bname)
    Vmarg = st_matrix(vname)
    m = rows(cells)
    if (p == 0 | any(X :>=.) | any(y :>=.) | any(beta :>=.))
        _jwidx_fail("Missing or empty nonlinear model design")
    if (sum(om)) {
        idx = selectindex(om':!=0)
        if (max(abs(bfull[,idx])) > 1e-12 | max(abs(vfull[idx,])) > 1e-12)
            _jwidx_fail("Omitted nonlinear model columns are not structural zeros")
    }
    eta = X*beta'
    if (max(abs(eta-st_data(.,xbvar))) > 1e-8*max((1,max(abs(eta)))))
        _jwidx_fail("Reconstructed nonlinear model columns disagree with the fitted index")
    mu = _jwidx_mean(eta,method)
    w = _jwidx_deriv(mu,method)
    if (any(mu :>=.) | any(w :>=.) | min(w)<=0)
        _jwidx_fail("The fitted nonlinear index has degenerate or overflowing derivatives")
    H = quadcross(X,w,X)
    Hinv = invsym(H)
    if (diag0cnt(Hinv) | any(Hinv:>=.))
        _jwidx_fail("The identified nonlinear design has a singular information matrix")
    err = max(abs(Hinv-MB))/max((1e-12,max(abs(MB))))
    if (err > 2e-5)
        _jwidx_fail("The reconstructed nonlinear Hessian disagrees with the original fit; refit jwdid")
    panel = panelsetup(id,1)
    N = rows(panel)
    if (N<=1) _jwidx_fail("At least two independent units are required")
    Q = panelsum(X:*(y-mu),panel)*Hinv
    scoresum = colsum(X:*(y-mu))'
    if (max(abs(scoresum):/(1:+colsum(abs(X:*(y-mu)))')) > 2e-5)
        _jwidx_fail("The nonlinear score is not centered at the stored fit; refit jwdid")
    V0 = quadcross(Q,Q)
    denom = sum(V0:^2)
    if (denom<=0 | missing(denom)) _jwidx_fail("Degenerate nonlinear score covariance")
    fac = sum(V0:*Vs)/denom
    err = max(abs(fac*V0-Vs))/max((1e-12,max(abs(Vs))))
    if (fac<=0 | missing(fac) | err>2e-5)
        _jwidx_fail("Unit-score covariance does not reproduce e(V) with a common finite-sample factor")

    stripe = st_matrixcolstripe(bbname)
    term = _jwidx_terms(stripe[active,2],gvar,tvar,xs)
    idx = selectindex(term[,1]:<.)
    if (!rows(idx)) _jwidx_fail("The fitted model contains no treatment-index coefficients")
    // Verify cohort support and the actual generated centered covariates.
    // This also detects stale centered variables before transporting scores.
    gc = uniqrows(sort(select(g,g:>0),1))
    for (j=1; j<=rows(gc); j++) {
        selectrow = selectindex(g:==gc[j])
        Z = (J(rows(selectrow),1,1), XS[selectrow,])
        if (rank(quadcross(Z,Z)) != q+1)
            _jwidx_fail("A treated cohort does not identify the full intercept-and-covariate index block")
        if (q) {
            meanx = mean(XS[selectrow,])
            for (k=1; k<=q; k++) {
                pick = selectindex((term[,1]:==gc[j]):*(term[,3]:==k))'
                if (!cols(pick)) _jwidx_fail("A treatment-index covariate component is absent")
                iscenter = term[pick[1],4]
                if (any(term[pick,4]:!=iscenter))
                    _jwidx_fail("Treatment-index centering differs within a cohort")
                if (iscenter) {
                    zcol = st_varindex("_x_"+xs[k])
                    if (missing(zcol)) _jwidx_fail("A generated centered covariate is unavailable")
                    tmp = st_data(selectrow,zcol)
                    if (max(abs(tmp-(XS[selectrow,k]:-meanx[k]))) >
                        1e-8*max((1,max(abs(XS[selectrow,k])))))
                        _jwidx_fail("A generated centered covariate differs from the original cohort mean centering")
                }
            }
        }
    }

    lead = selectindex(cells[,2]:<(cells[,1]:-gap*anti:-1e-8*gap))
    nlead = rows(lead)
    K = nlead*(q+1)
    L = J(K,p,0)
    HM = J(K,2,.)
    IT = J(K,2,.)
    A = J(m,p,0)
    D = J(N,m,0)
    X0 = X
    X0[,idx] = J(n,rows(idx),0)
    eta0 = X0*beta'
    mu0 = _jwidx_mean(eta0,method)
    w0 = _jwidx_deriv(mu0,method)
    if (any(mu0:>=.) | any(w0:>=.))
        _jwidx_fail("Counterfactual nonlinear predictions overflow")
    r = 0
    for (j=1; j<=m; j++) {
        gval = cells[j,1]
        tval = cells[j,2]
        selectrow = selectindex((g:==gval):*(t:==tval))
        ng = rows(selectrow)
        if (!ng) _jwidx_fail("A source cohort-time margin has no observations")
        if (abs(tval-(gval-gap*anti))<1e-8*gap) {
            if (abs(b[j])>1e-10) _jwidx_fail("The normalized baseline ATT is nonzero")
            continue
        }
        coeff = J(q+1,1,.)
        for (k=0; k<=q; k++) {
            src = selectindex((term[,1]:==gval):*(term[,2]:==tval):*(term[,3]:==k))
            if (rows(src)!=1)
                _jwidx_fail("A cohort-time index block is missing, omitted, or duplicated; refit an identified full block")
            coeff[k+1] = src[1]
        }
        f = mu[selectrow]-mu0[selectrow]
        val = mean(f)
        if (abs(val-b[j])>2e-8*max((1,abs(b[j]))))
            _jwidx_fail("Analytic outcome contrasts disagree with jwdid's margins")
        A[j,] = mean(X[selectrow,]:*w[selectrow]-X0[selectrow,]:*w0[selectrow])
        if (vce == "unconditional") {
            // One row per cohort member at this period. In sampling-contribution
            // units this is D_g/N_g * (f_i-d), not D_g/N * (f_i-d).
            D[id[selectrow],j] = (f:-val)/ng
        }
        if (tval < gval-gap*anti-1e-8*gap) {
            meanx = J(1,q,0)
            if (q) meanx = mean(XS[selectrow,])
            L[r+1,coeff[1]] = 1
            for (k=1; k<=q; k++) {
                if (term[coeff[k+1],4]) L[r+1,coeff[k+1]] = -meanx[k]
                L[r+k+1,coeff[k+1]] = 1
            }
            HM[|r+1,1\r+q+1,2|] = J(q+1,1,1)*(gval,tval)
            IT[|r+1,1\r+q+1,2|] = (active[coeff]',(0::q))
            r = r+q+1
        }
    }
    if (r!=K) _jwidx_fail("Internal index-block dimension mismatch")
    QD = Q*A'
    margerr = max(abs(fac*quadcross(QD,QD)-Vmarg))/max((1e-12,max(abs(Vmarg))))
    if (margerr>5e-5)
        _jwidx_fail("Analytic conditional ATT covariance disagrees with jwdid's numerical margins Jacobian")
    QD = QD+D
    if (K) {
        h = beta*L'
        QH = Q*L'
        joint = (QD,QH)
        VV = fac*quadcross(joint,joint)
        st_matrix(hname,h)
        st_matrix(vhname,VV[|m+1,m+1\m+K,m+K|])
        st_matrix(vdhname,VV[|1,m+1\m,m+K|])
        st_matrix(hmname,HM)
        st_matrix(itname,IT)
        st_matrix(vname,VV[|1,1\m,m|])
    }
    else st_matrix(vname,fac*quadcross(QD,QD))
    st_numscalar(nkname,K)
    st_numscalar(facname,fac)
    st_numscalar(errname,err)
    st_numscalar(merrname,margerr)
}
end
