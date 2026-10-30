*! jwselect 0.7.0 30oct2026: read-only linear/nonlinear jwdid adapter
program define _jwselect_jwdid, rclass
    version 15.1
    syntax [, VCE(string)]

    if `"`vce'"' == "" local vce unconditional
    local vce = lower(`"`vce'"')
    if !inlist(`"`vce'"', "conditional", "unconditional") {
        di as err "vce() must be conditional or unconditional"
        exit 198
    }
    if `"`e(cmd)'"' != "jwdid" {
        di as err "jwselect requires active jwdid results or the matrix interface"
        exit 301
    }

    tempname held B V C T
    quietly estimates store `held'
    preserve
    capture noisily _jwselect_jwdid_work, vce(`vce')
    local rc = _rc
    if !`rc' {
        matrix `B' = r(bgt)
        matrix `V' = r(Vgt)
        matrix `C' = r(cells)
        matrix `T' = r(upstream_table)
        foreach key in N_units N_obs N_clusters N_covariates gap anticipation df_r {
            local `key' = r(`key')
        }
        foreach key in ivar tvar gvar sourcecmd sourcevce upstream_stat model method restriction_scale {
            local `key' `"`r(`key')'"'
        }
    }
    restore
    quietly estimates restore `held'
    quietly estimates drop `held'
    if `rc' exit `rc'

    return matrix bgt = `B'
    return matrix Vgt = `V'
    return matrix cells = `C'
    return matrix upstream_table = `T'
    foreach key in N_units N_obs N_clusters N_covariates gap anticipation df_r {
        return scalar `key' = ``key''
    }
    foreach key in ivar tvar gvar sourcecmd sourcevce upstream_stat model method restriction_scale {
        return local `key' `"``key''"'
    }
    return local vce `"`vce'"'
end

program define _jwselect_jwdid_work, rclass
    version 15.1
    syntax , VCE(string)

    local backend `"`e(cmd2)'"'
    if !inlist(`"`backend'"', "", "regress", "logit", "poisson") {
        di as err "The jwdid adapter supports regress, default reghdfe, logit, and poisson"
        di as err "Other methods require a separately validated matrix interface"
        exit 198
    }
    local model linear
    local method `backend'
    if `"`method'"' == "" local method reghdfe
    local restriction_scale outcome_contrast
    if inlist(`"`backend'"', "logit", "poisson") {
        local model nonlinear
        local restriction_scale outcome_placebo
        if missing(e(converged)) | e(converged) != 1 {
            di as err "The nonlinear jwdid fit did not converge"
            exit 498
        }
        if `"`e(index)'"' != "" {
            di as err "The nonlinear adapter requires outcome-scale margins; refit without index"
            exit 198
        }
    }
    if `"`e(cmdopt)'"' != "" {
        di as err "Custom method() options are outside this adapter's validated contract"
        exit 198
    }
    if `"`e(hettype)'"' != "timecohort" | `"`e(type)'"' != "never" {
        di as err "The jwdid adapter requires hettype(timecohort) and never"
        exit 198
    }
    if `"`vce'"' == "unconditional" & `"`backend'"' == "" {
        di as err "vce(unconditional) requires jwdid ..., method(regress, logit, or poisson) never"
        di as err "Use vce(conditional) for an explicitly conditional analysis of a reghdfe fit"
        exit 198
    }

    local ivar `"`e(ivar)'"'
    local tvar `"`e(tvar)'"'
    local gvar `"`e(gvar)'"'
    local sourcecmd `"`e(cmdline)'"'
    local sourcevce `"`e(vce)'"'
    local sourceN = e(N)
    local gap = e(gap)
    local anti = e(anticipation)
    local ncl = .
    local df = .
    capture local ncl = e(N_clust)
    capture local df = e(df_r)
    if `"`ivar'"' == "" | `"`tvar'"' == "" | `"`gvar'"' == "" {
        di as err "The jwdid adapter requires ivar(), tvar(), and gvar() metadata"
        exit 198
    }
    confirm variable `ivar'
    confirm numeric variable `tvar' `gvar'
    confirm numeric variable __tr__
    if missing(`gap') | `gap' <= 0 | missing(`anti') | `anti' <= 0 | `anti' != floor(`anti') {
        di as err "Missing or invalid jwdid period-gap / anticipation metadata"
        exit 498
    }

    // Parse the original command, not the expanded regression formula.
    gettoken cmd commandargs : sourcecmd
    if lower(`"`cmd'"') != "jwdid" {
        di as err "Cannot validate this fit: e(cmdline) is not an original jwdid command"
        exit 498
    }
    quietly _jwselect_jwdid_parse `commandargs'
    local xs `"`r(covariates)'"'
    local y `"`r(outcome)'"'
    local ncov = r(N_covariates)
    local cluster `"`r(cluster)'"'
    if `"`cluster'"' != "" & `"`cluster'"' != `"`ivar'"' {
        di as err "The jwdid adapter requires clustering at ivar(), its independent sampling unit"
        di as err "Use the matrix interface with an explicit independent-unit count for other designs"
        exit 198
    }

    tempvar sample idtag timetag ni
    quietly generate byte `sample' = e(sample)
    quietly count if `sample'
    if r(N) == 0 | r(N) != `sourceN' {
        di as err "The original estimation sample is unavailable or has changed"
        exit 498
    }
    // The outer preserve restores all observations, labels, variables, and sort order.
    // Restricting here also prevents out-of-sample cells from entering margins.
    quietly keep if `sample'
    if `"`backend'"' == "logit" {
        capture assert inlist(`y', 0, 1)
        if _rc {
            di as err "The logit adapter requires a binary outcome coded 0/1"
            exit 198
        }
    }
    capture isid `ivar' `tvar'
    if _rc {
        di as err "The jwdid adapter requires one observation per individual and period"
        exit 459
    }
    quietly assert !missing(`gvar') & `gvar' >= 0
    quietly bysort `ivar' (`tvar'): assert `gvar' == `gvar'[1]
    quietly count if `gvar' == 0
    if r(N) == 0 {
        di as err "An observed never-treated cohort coded gvar=0 is required"
        di as err "jwdid's last-cohort fallback is outside the supported adapter contract"
        exit 198
    }
    foreach x of local xs {
        capture by `ivar': assert `x' == `x'[1]
        if _rc {
            di as err "Covariate `x' varies within individual; this adapter requires time-invariant covariates"
            exit 198
        }
    }
    quietly egen byte `idtag' = tag(`ivar')
    quietly count if `idtag'
    local nu = r(N)
    if missing(`ncl') | `ncl' != `nu' {
        di as err "The recorded cluster count does not match the number of individuals"
        exit 498
    }
    quietly egen byte `timetag' = tag(`tvar')
    quietly count if `timetag'
    local nt = r(N)
    quietly bysort `ivar': generate long `ni' = _N
    capture assert `ni' == `nt'
    if _rc {
        di as err "The jwdid adapter requires a balanced panel with common time support"
        exit 459
    }
    quietly levelsof `tvar', local(times)
    local previous
    foreach tt of local times {
        if `"`previous'"' != "" {
            if abs((`tt' - `previous') - `gap') > 1e-8 * max(1, abs(`gap')) {
                di as err "Observed time support is not regularly spaced by e(gap)"
                exit 459
            }
        }
        local previous `tt'
    }

    tempname B V C T miss asym
    return clear
    if `"`vce'"' == "unconditional" {
        capture quietly estat attgt, vce(unconditional)
    }
    else {
        capture quietly estat attgt
    }
    local marginrc = _rc
    if `marginrc' {
        di as err "jwdid cohort-time margins extraction failed (return code `marginrc')"
        di as err "Restore the original jwdid fit and run estat attgt directly to inspect the upstream error"
        exit `marginrc'
    }
    // jwdid_estat captures its own errors. Do not trust its zero exit code alone.
    if `"`r(agg)'"' != "attgt" | `"`r(ecmd)'"' != "jwdid_estat" {
        di as err "jwdid did not return fresh cohort-time margins"
        di as err "Restore the original jwdid fit and run estat attgt directly to inspect the upstream error"
        exit 498
    }
    capture matrix `T' = r(table)
    if _rc {
        di as err "jwdid did not return its cohort-time inference table"
        exit 498
    }
    local tablerows : rownames `T'
    local upstream_stat : word 3 of `tablerows'
    capture matrix `B' = r(b)
    if _rc {
        di as err "jwdid did not return the cohort-time estimate vector"
        exit 498
    }
    capture matrix `V' = r(V)
    if _rc {
        di as err "jwdid did not return the full joint cohort-time covariance matrix"
        exit 498
    }
    if rowsof(`T') < 6 | colsof(`T') != colsof(`B') {
        di as err "Invalid dimensions in jwdid's cohort-time inference table"
        exit 503
    }
    if rowsof(`B') != 1 | colsof(`B') == 0 | ///
        rowsof(`V') != colsof(`B') | colsof(`V') != colsof(`B') {
        di as err "Invalid dimensions in jwdid's cohort-time returns"
        exit 503
    }
    mata: st_numscalar("`miss'", sum(st_matrix("`B'"):>=.) + sum(st_matrix("`V'"):>=.))
    if scalar(`miss') {
        di as err "Cohort-time margins contain missing estimates or covariance entries"
        exit 498
    }
    mata: st_numscalar("`asym'", max(abs(st_matrix("`V'") - st_matrix("`V'")')) / max((1,max(abs(st_matrix("`V'"))))))
    if scalar(`asym') > 1e-8 {
        di as err "jwdid returned a nonsymmetric cohort-time covariance matrix"
        exit 498
    }
    confirm numeric variable __attgt__
    local m = colsof(`B')
    matrix `C' = J(`m', 2, .)
    local names : colnames `B'
    local keep
    forvalues j = 1/`m' {
        local cname : word `j' of `names'
        if !regexm(`"`cname'"', "^([0-9]+)[bon]*[.]__attgt__$") {
            di as err "Unrecognized jwdid cohort-time margin name: `cname'"
            exit 498
        }
        local code = real(regexs(1))
        quietly count if __attgt__ == `code'
        if r(N) == 0 {
            di as err "Cannot map margin `cname' to the original estimation sample"
            exit 498
        }
        quietly summarize `gvar' if __attgt__ == `code', meanonly
        if r(min) != r(max) | missing(r(min)) {
            di as err "A cohort-time margin maps to more than one cohort"
            exit 498
        }
        matrix `C'[`j',1] = r(min)
        if r(min) > 0 local keep `keep' `j'
        quietly summarize `tvar' if __attgt__ == `code', meanonly
        if r(min) != r(max) | missing(r(min)) {
            di as err "A cohort-time margin maps to more than one period"
            exit 498
        }
        matrix `C'[`j',2] = r(min)
    }
    if `"`keep'"' == "" {
        di as err "No treated cohort-time margins were returned"
        exit 498
    }
    mata: st_matrix("`T'", st_matrix("`T'")[1..6,strtoreal(tokens(st_local("keep"))) ]')
    // Retain structural zero reference-period rows. Their singularity is expected.
    mata: st_matrix("`B'", st_matrix("`B'")[,strtoreal(tokens(st_local("keep")))])
    mata: st_matrix("`V'", st_matrix("`V'")[strtoreal(tokens(st_local("keep"))),strtoreal(tokens(st_local("keep")))])
    mata: st_matrix("`C'", st_matrix("`C'")[strtoreal(tokens(st_local("keep"))),])
    // Nonlinear effects are outcome-scale predictive contrasts. The core may
    // use their zero-placebo restrictions, but must not rebase them by
    // subtracting a nonzero pretreatment margin. Main-command guards enforce
    // the original cohort baseline and anticipation. Verify its normalization.
    // margins can leave roundoff in this identically-zero contrast because
    // it evaluates numerical derivatives. Only after the 1e-10 checks below,
    // canonicalize this known normalization to exact zero. Materially nonzero
    // estimates/covariances are errors, never silently removed.
    if `"`model'"' == "nonlinear" {
        tempname invalidbase
        quietly levelsof `gvar' if `gvar'>0, local(groups)
        foreach g of local groups {
            local found 0
            forvalues j=1/`=rowsof(`C')' {
                if `C'[`j',1]==`g' & abs(`C'[`j',2]-(`g'-`gap'*`anti'))<1e-8*`gap' {
                    local ++found
                    mata: st_numscalar("`invalidbase'", max(abs((st_matrix("`V'")[`j',.],st_matrix("`V'")[.,`j']'))))
                    if abs(`B'[1,`j'])>1e-10 | scalar(`invalidbase')>1e-10 {
                        di as err "Nonlinear cohort baseline is not a structural-zero margin"
                        exit 498
                    }
                    matrix `B'[1,`j']=0
                    matrix `V'[`j',1]=J(1,colsof(`V'),0)
                    matrix `V'[1,`j']=J(rowsof(`V'),1,0)
                    // Structural zeros have no meaningful test statistic,
                    // p-value, or interval in the upstream comparison table.
                    matrix `T'[`j',3]=J(1,4,.)
                }
            }
            if `found'!=1 {
                di as err "Nonlinear cohort baseline is unavailable or duplicated"
                exit 498
            }
        }
    }
    // r(table) can leave normalized rows missing. Use the validated estimate
    // vector and covariance for the displayed estimates and standard errors.
    // Retain upstream test statistics, p-values and interval endpoints as given.
    mata: st_matrix("`T'", (st_matrix("`B'")', sqrt(diagonal(st_matrix("`V'"))), st_matrix("`T'")[.,3..6]))
    matrix colnames `C' = cohort time
    matrix `T' = `C', `T'
    matrix colnames `T' = cohort time estimate se statistic pvalue ll ul

    return matrix bgt = `B'
    return matrix Vgt = `V'
    return matrix cells = `C'
    return matrix upstream_table = `T'
    return local upstream_stat `upstream_stat'
    return scalar N_units = `nu'
    return scalar N_obs = `sourceN'
    return scalar N_clusters = `ncl'
    return scalar N_covariates = `ncov'
    return scalar gap = `gap'
    return scalar anticipation = `anti'
    return scalar df_r = `df'
    return local ivar `"`ivar'"'
    return local tvar `"`tvar'"'
    return local gvar `"`gvar'"'
    return local sourcecmd `"`sourcecmd'"'
    return local sourcevce `"`sourcevce'"'
    return local model `model'
    return local method `method'
    return local restriction_scale `restriction_scale'
end

program define _jwselect_jwdid_parse, rclass
    version 15.1
    syntax varlist(numeric fv ts) [if] [in] [pw iw aw], ///
        [Ivar(varname) cluster(varlist) Tvar(varname) time(varname) ///
        fevar(varlist fv ts) Gvar(varname) trtvar(varname) trgvar(varname) ///
        never group method(string asis) corr hettype(string asis) ///
        exovar(varlist fv ts) exogvar(varlist fv ts) xtvar(varlist fv ts) ///
        xgvar(varlist fv ts) xattvar(varlist fv ts) xasis cre ///
        ANTIcipation(numlist max=1 >0) index *]
    if `"`weight'"' != "" {
        di as err "The initial jwdid adapter supports unweighted fits; use bmat()/vmat()/cells() for other validated inputs"
        exit 198
    }
    if `"`corr'`cre'`trtvar'`trgvar'`fevar'`exovar'`exogvar'`xtvar'`xgvar'`xattvar'"' != "" {
        di as err "This jwdid specification uses a custom mode outside the adapter's validated contract"
        di as err "Unsupported modes: corr, cre, trtvar, trgvar, fevar, exovar, exogvar, xtvar, xgvar, xattvar"
        di as err "Use the matrix interface after validating the required cohort-time estimand and covariance"
        exit 198
    }
    if `"`options'"' != "" {
        di as err "Extra regression options are outside the initial jwdid adapter contract: `options'"
        exit 198
    }
    gettoken y xs : varlist
    capture confirm numeric variable `y'
    if _rc {
        di as err "This adapter requires a plain numeric outcome; generate transformed outcomes explicitly before fitting jwdid"
        exit 198
    }
    foreach x of local xs {
        capture confirm numeric variable `x'
        if _rc {
            di as err "This adapter accepts plain numeric covariates; expand factor/TS expressions explicitly before fitting jwdid"
            exit 198
        }
    }
    return local covariates `"`xs'"'
    return local outcome `"`y'"'
    return local cluster `"`cluster'"'
    return scalar N_covariates = `: word count `xs''
end
