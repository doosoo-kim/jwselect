*! jwselect 0.7.0 30oct2026: read-only CS-DR adapter for csdid 2.0.1+
program define _jwselect_csdid, rclass
    version 16.1
    syntax [, VCE(string)]
    if "`vce'" == "" local vce supplied
    if lower("`vce'") != "supplied" {
        di as err "The csdid adapter uses vce(supplied): csdid's analytical influence-function covariance"
        exit 198
    }
    if "`e(cmd)'" != "csdid" {
        di as err "The CS-DR adapter requires active csdid results"
        exit 301
    }

    tempname held B V C T
    quietly estimates store `held'
    preserve
    capture noisily _jwselect_csdid_work
    local rc = _rc
    if !`rc' {
        matrix `B' = r(bgt)
        matrix `V' = r(Vgt)
        matrix `C' = r(cells)
        matrix `T' = r(upstream_table)
        foreach key in N_units N_obs N_clusters N_covariates gap anticipation df_r pscoretrim {
            local `key' = r(`key')
        }
        foreach key in ivar tvar gvar sourcecmd sourcevce sourceversion sourceestimator {
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
    foreach key in N_units N_obs N_clusters N_covariates gap anticipation df_r pscoretrim {
        return scalar `key' = ``key''
    }
    foreach key in ivar tvar gvar sourcecmd sourcevce sourceversion sourceestimator {
        return local `key' `"``key''"'
    }
    return local vce supplied
    return local model csdr
    return local method dr
    return local restriction_scale outcome_placebo
    return local upstream_stat z
end

program define _jwselect_csdid_work, rclass
    version 16.1
    // Only the documented 2.x ATT(g,t) interface is supported. In particular,
    // a legacy e(V_attgt) and a 2.x bootstrap e(V) are not interchangeable.
    local sourceversion `"`e(version)'"'
    local compatible 0
    if regexm("`sourceversion'", "^2[.]([0-9]+)[.]([0-9]+)$") {
        local minor = real(regexs(1))
        local patch = real(regexs(2))
        local compatible = (`minor' > 0 | `patch' >= 1)
    }
    if !`compatible' {
        di as err "jwselect's CS-DR adapter requires csdid 2.0.1 or a compatible 2.x release"
        di as err "Legacy SSC csdid 1.x results are not supported; use the documented 2.x installation"
        exit 198
    }
    if "`e(method)'" != "dr" {
        di as err "The CS-DR adapter requires csdid ..., method(dr)"
        exit 198
    }
    if "`e(control_group)'" != "nevertreated" | "`e(base_period)'" != "universal" {
        di as err "The initial CS-DR adapter requires nevertreated base_period(universal)"
        di as err "Refit csdid with those options; jwselect does not change the comparison group or rebase DR estimates"
        exit 198
    }
    if e(bstrap) != 0 | !inlist("`e(vce)'", "analytical", "cluster") {
        di as err "The CS-DR adapter requires csdid's analytical covariance; refit with analytical"
        exit 198
    }
    if "`e(panel_mode)'" != "panel" {
        di as err "The initial CS-DR adapter requires a balanced panel (ivar() and bal(full))"
        exit 198
    }
    if `"`e(wtype)'`e(wexp)'`e(weightvar)'"' != "" {
        di as err "The initial CS-DR adapter supports unweighted fits only"
        exit 198
    }
    local ivar `"`e(idvar)'"'
    local tvar `"`e(timevar)'"'
    local gvar `"`e(gvar)'"'
    local clvar `"`e(clustervar)'"'
    local sourcecmd `"`e(cmdline)'"'
    local sourcevce `"`e(vce)'"'
    local sourceN = e(N)
    local nu = e(N_units)
    local nt = e(N_time)
    local ng = e(N_groups)
    local anti = e(anticipation)
    local trim = e(pscoretrim)
    local lev = e(level)
    if "`ivar'" == "" | "`tvar'" == "" | "`gvar'" == "" {
        di as err "csdid panel, time, and treatment-cohort metadata are required"
        exit 498
    }
    if "`clvar'" != "" & "`clvar'" != "`ivar'" {
        di as err "The initial CS-DR adapter requires independent panel units; cluster() must be ivar() or omitted"
        exit 198
    }
    confirm numeric variable `ivar' `tvar' `gvar'
    if missing(`nu', `sourceN', `nt', `ng', `anti', `lev', `trim') | ///
        `nu' <= 1 | `nt' < 2 | `ng' < 1 | `anti' < 0 | `anti' != floor(`anti') | ///
        `sourceN' != `nu' * `nt' | `lev' <= 0 | `lev' >= 100 {
        di as err "Invalid csdid sample, anticipation, or inference metadata"
        exit 498
    }
    if "`clvar'" != "" {
        if e(N_clusters) != `nu' {
            di as err "csdid's cluster count does not match its panel-unit count"
            exit 498
        }
    }
    // csdid signs outcome, covariates, identifiers, and its estimation sample.
    // Require that signature: equal observation counts alone do not detect edits.
    if `"`e(datasignaturevars)'"' == "" {
        di as err "The original signed csdid estimation sample is unavailable; restore it and refit csdid"
        exit 498
    }
    capture checkestimationsample
    if _rc {
        di as err "The data have changed since csdid was fitted; restore the estimation data or refit csdid"
        exit 459
    }
    gettoken cmd commandargs : sourcecmd
    if lower("`cmd'") != "csdid" {
        di as err "Cannot validate the original csdid command"
        exit 498
    }
    quietly _jwselect_csdid_parse `commandargs'
    local xs `"`r(covariates)'"'
    local ncov = r(N_covariates)

    tempvar sample idtag timetag ni gtag
    quietly generate byte `sample' = e(sample)
    quietly count if `sample'
    if r(N) != `sourceN' {
        di as err "The csdid estimation sample is unavailable or its size has changed"
        exit 498
    }
    quietly keep if `sample'
    capture isid `ivar' `tvar'
    if _rc {
        di as err "The CS-DR adapter requires one observation per individual and period"
        exit 459
    }
    capture assert !missing(`ivar', `gvar', `tvar') & `gvar' >= 0 & ///
        `gvar' == floor(`gvar') & `tvar' > 0 & `tvar' == floor(`tvar')
    if _rc {
        di as err "The initial CS-DR adapter requires positive integer calendar times and nonnegative integer cohorts"
        exit 198
    }
    capture bysort `ivar' (`tvar'): assert `gvar' == `gvar'[1]
    if _rc {
        di as err "Treatment cohort must be constant within panel unit"
        exit 459
    }
    foreach x of local xs {
        capture by `ivar': assert `x' == `x'[1]
        if _rc {
            di as err "Covariate `x' varies within panel unit; the initial CS-DR adapter requires time-invariant covariates"
            exit 198
        }
    }
    quietly count if `gvar' == 0
    if r(N) == 0 {
        di as err "An observed never-treated group coded gvar=0 is required"
        exit 198
    }
    quietly egen byte `idtag' = tag(`ivar')
    quietly count if `idtag'
    if r(N) != `nu' {
        di as err "The number of observed panel units differs from csdid's metadata"
        exit 498
    }
    quietly egen byte `timetag' = tag(`tvar')
    quietly count if `timetag'
    if r(N) != `nt' {
        di as err "Observed time support differs from csdid's metadata"
        exit 498
    }
    quietly bysort `ivar': generate long `ni' = _N
    capture assert `ni' == `nt'
    if _rc {
        di as err "The CS-DR adapter requires a balanced panel with common time support"
        exit 459
    }
    quietly levelsof `tvar', local(times)
    local previous
    foreach tt of local times {
        if "`previous'" != "" {
            if `tt' - `previous' != 1 {
                di as err "The initial CS-DR adapter requires unit-spaced calendar times"
                exit 198
            }
        }
        local previous `tt'
    }
    quietly levelsof `gvar' if `gvar' > 0, local(groups)
    if `: word count `groups'' != `ng' {
        di as err "Observed treatment cohorts differ from csdid's metadata"
        exit 498
    }

    tempname A B0 V0 B V C T L miss asym eig
    foreach mat in attgt b V {
        capture confirm matrix e(`mat')
        if _rc {
            di as err "csdid did not return the required e(`mat') matrix"
            exit 498
        }
    }
    matrix `A' = e(attgt)
    matrix `B0' = e(b)
    matrix `V0' = e(V)
    local m = rowsof(`A')
    local k = colsof(`B0')
    if colsof(`A') != 10 | `m' != `ng' * `nt' | rowsof(`B0') != 1 | ///
        `k' != `m' - `ng' | rowsof(`V0') != `k' | colsof(`V0') != `k' {
        di as err "Expected a complete unaggregated csdid ATT(g,t) grid and its joint covariance"
        di as err "Restore the original csdid estimates before aggregation posting; failed or missing cells are not supported"
        exit 498
    }
    local acols : colnames `A'
    if "`acols'" != "group time event_time att se n_treat_t n_treat_pre n_control_t n_control_pre base_time" {
        di as err "Unrecognized csdid e(attgt) column layout"
        exit 498
    }
    local bnames : colnames `B0'
    local vnames : colnames `V0'
    local vrows : rownames `V0'
    if "`bnames'" != "`vnames'" | "`bnames'" != "`vrows'" {
        di as err "csdid's coefficient and covariance names are not aligned"
        exit 498
    }
    // csdid intentionally prints missing se for its structural-zero reference
    // row. Only those se entries may be missing; validate nonbase se below.
    mata: st_numscalar("`miss'", hasmissing(st_matrix("`A'")[.,(1..4,6..10)]) | hasmissing(st_matrix("`B0'")) | hasmissing(st_matrix("`V0'")))
    if scalar(`miss') {
        di as err "csdid returned missing estimates, standard errors, or covariance entries"
        exit 498
    }
    mata: st_numscalar("`asym'", max(abs(st_matrix("`V0'") - st_matrix("`V0'")')) / max((1,max(abs(st_matrix("`V0'"))))))
    if scalar(`asym') > 1e-10 {
        di as err "csdid returned a nonsymmetric joint covariance"
        exit 498
    }
    mata: st_numscalar("`eig'", min(symeigenvalues(st_matrix("`V0'"))) / max((1,max(abs(st_matrix("`V0'"))))))
    if scalar(`eig') < -1e-10 {
        di as err "csdid's joint covariance is not positive semidefinite"
        exit 498
    }
    matrix `L' = J(`m', `k', 0)
    matrix `C' = `A'[1...,1..2]
    local used
    local normcount 0
    forvalues i = 1/`m' {
        local gg = `A'[`i',1]
        local tt = `A'[`i',2]
        local bb = `A'[`i',10]
        local gmatch : list gg in groups
        local tmatch : list tt in times
        if !`gmatch' | !`tmatch' | `bb' != `gg' - `anti' - 1 | ///
            `A'[`i',3] != `tt' - `gg' {
            di as err "csdid cell metadata do not match the declared universal baseline and sample"
            exit 498
        }
        if `i' > 1 {
            forvalues j = 1/`=`i'-1' {
                if `C'[`j',1] == `gg' & `C'[`j',2] == `tt' {
                    di as err "csdid returned duplicate cohort-time cells"
                    exit 498
                }
            }
        }
        if `tt' == `bb' {
            if `A'[`i',4] != 0 | (!missing(`A'[`i',5]) & `A'[`i',5] != 0) {
                di as err "csdid's universal-base normalization must have exactly zero estimate and variance"
                exit 498
            }
            local ++normcount
        }
        else {
            // Integer calendar support makes this documented name unambiguous.
            // Search by name so order changes in e(b) or e(attgt) stay harmless.
            local cname g`gg'___`tt'_`bb'
            local pos : list posof "`cname'" in bnames
            if !`pos' {
                di as err "Cannot map csdid cell `cname' to e(b); restore the unaggregated ATT(g,t) results"
                exit 498
            }
            local duplicate : list pos in used
            if `duplicate' {
                di as err "csdid's coefficient names do not identify unique cohort-time cells"
                exit 498
            }
            local used `used' `pos'
            if abs(`B0'[1,`pos'] - `A'[`i',4]) > 1e-10 * max(1,abs(`A'[`i',4])) | ///
                missing(`A'[`i',5]) | `A'[`i',5] <= 0 | ///
                abs(`V0'[`pos',`pos'] - `A'[`i',5]^2) > 1e-10 * max(1,`A'[`i',5]^2) {
                di as err "csdid's posted coefficients/covariance do not match its ATT(g,t) table"
                exit 498
            }
            matrix `L'[`i',`pos'] = 1
        }
    }
    if `normcount' != `ng' | `: word count `used'' != `k' {
        di as err "Each cohort must have one structural-zero baseline and every posted coefficient must map once"
        exit 498
    }
    matrix `B' = `B0' * `L''
    matrix `V' = `L' * `V0' * `L''
    matrix `T' = J(`m', 8, .)
    matrix `T'[1,1] = `C'
    local critical = invnormal(1-(100-`lev')/200)
    forvalues i = 1/`m' {
        matrix `T'[`i',3] = `B'[1,`i']
        matrix `T'[`i',4] = sqrt(`V'[`i',`i'])
        if `T'[`i',4] > 0 {
            matrix `T'[`i',5] = `T'[`i',3] / `T'[`i',4]
            matrix `T'[`i',6] = 2*normal(-abs(`T'[`i',5]))
            matrix `T'[`i',7] = `T'[`i',3] - `critical' * `T'[`i',4]
            matrix `T'[`i',8] = `T'[`i',3] + `critical' * `T'[`i',4]
        }
    }
    matrix colnames `C' = cohort time
    matrix colnames `T' = cohort time estimate se statistic pvalue ll ul
    return matrix bgt = `B'
    return matrix Vgt = `V'
    return matrix cells = `C'
    return matrix upstream_table = `T'
    return scalar N_units = `nu'
    return scalar N_obs = `sourceN'
    return scalar N_clusters = `nu'
    return scalar N_covariates = `ncov'
    return scalar gap = 1
    return scalar anticipation = `anti'
    return scalar df_r = .
    return scalar pscoretrim = `trim'
    return local ivar "`ivar'"
    return local tvar "`tvar'"
    return local gvar "`gvar'"
    return local sourcecmd `"`sourcecmd'"'
    return local sourcevce "`sourcevce'"
    return local sourceversion "`sourceversion'"
    return local sourceestimator csdid
end

program define _jwselect_csdid_parse, rclass
    version 16.1
    syntax varlist(numeric fv) [if] [in] [iw], [*]
    if "`weight'" != "" {
        di as err "The initial CS-DR adapter supports unweighted fits only"
        exit 198
    }
    gettoken y xs : varlist
    capture confirm numeric variable `y'
    if _rc {
        di as err "Use a plain numeric outcome in the csdid fit"
        exit 198
    }
    foreach x of local xs {
        capture confirm numeric variable `x'
        if _rc {
            di as err "The initial CS-DR adapter accepts plain numeric covariates; expand factor expressions before fitting"
            exit 198
        }
    }
    return local covariates "`xs'"
    return scalar N_covariates = `: word count `xs''
end
