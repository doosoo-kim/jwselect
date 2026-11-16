*! jwselect 0.9.1 16nov2026: population interaction-weighted SA source
program define jwselect_sa_iw, eclass
    version 16.1
    if `"`0'"'=="" {
        if "`e(cmd)'"!="jwselect_sa_iw" error 301
        if `"`0'"'!="" exit 198
        _jwselect_sa_iw_display
        exit
    }
    local original0 `"`0'"'
    syntax, POSTWindow(numlist integer min=2 max=2) ///
        [COHORTS(numlist integer sort) ANTicipation(integer 0) BASELINE(string) Level(cilevel)]
    if "`baseline'"=="" local baseline cohort
    if !inlist("`baseline'","cohort","common") | `anticipation'<0 {
        di as err "baseline() must be cohort or common and anticipation() nonnegative"
        exit 198
    }
    local lo : word 1 of `postwindow'
    local hi : word 2 of `postwindow'
    if `lo'<0 | `hi'<`lo' {
        di as err "postwindow() requires increasing nonnegative event-time endpoints"
        exit 198
    }
    // Validates the full source, signed data, support and upstream map.
    _jwselect_sa, vce(supplied)
    local ivar `r(ivar)'
    local tvar `r(tvar)'
    local gvar `r(gvar)'
    local y `e(sa_outcome)'
    local nu=e(sa_N_units)
    local nt=e(sa_N_time)
    local sourcecmd `"`e(cmdline)'"'
    tempname previous D H V HM W TW BASES TARGETMAP K sharev fixedv B VD
    tempvar touse
    quietly generate byte `touse'=e(sample)
    quietly summarize `tvar' if `touse', meanonly
    local first=r(min)
    estimates store `previous'
    capture mata: assert(jw_sa_iw_version()=="0.9.1")
    if _rc {
        quietly findfile _jwselect_sa_iw_mata.do
        quietly do `"`r(fn)'"'
    }
    preserve
    quietly sort `ivar' `tvar'
    capture noisily mata: jw_sa_iw_inputs("`y'","`gvar'","`touse'", ///
        `nt',`first',strtoreal(tokens("`postwindow'")), ///
        strtoreal(tokens("`cohorts'"))',`anticipation',"`baseline'", ///
        "`D'","`H'","`V'","`HM'","`W'","`TW'","`BASES'","`TARGETMAP'", ///
        "`K'","`sharev'","`fixedv'")
    local rc=_rc
    restore
    if `rc' {
        quietly estimates restore `previous'
        quietly estimates drop `previous'
        exit `rc'
    }
    matrix `B'=`D'
    matrix colnames `B'=population_iw
    matrix `VD'=`V'[1,1]
    matrix rownames `VD'=population_iw
    matrix colnames `VD'=population_iw
    local nobs=`nu'*`nt'
    ereturn post `B' `VD', esample(`touse') obs(`nobs')
    ereturn matrix jw_d=`D'
    ereturn matrix jw_V=`V'
    if scalar(`K')>0 {
        ereturn matrix jw_h=`H'
        ereturn matrix jw_hmap=`HM'
    }
    ereturn matrix jw_bases=`BASES'
    ereturn matrix jw_targetmap=`TARGETMAP'
    matrix colnames `W'=cohort N share weight baseline_effect
    matrix colnames `TW'=cohort time weight
    ereturn matrix iw_cohorts=`W'
    ereturn matrix iw_targetweights=`TW'
    ereturn scalar iw_share_variance=scalar(`sharev')
    ereturn scalar iw_fixed_variance=scalar(`fixedv')
    ereturn scalar jw_nunits=`nu'
    ereturn scalar jw_K=scalar(`K')
    ereturn scalar N_units=`nu'
    ereturn scalar N_clusters=`nu'
    ereturn scalar N_obs=`nu'*`nt'
    ereturn scalar N_covariates=0
    ereturn scalar gap=1
    ereturn scalar anticipation=`anticipation'
    ereturn scalar level=`level'
    ereturn local baseline "`baseline'"
    ereturn local postwindow "`postwindow'"
    ereturn local ivar "`ivar'"
    ereturn local tvar "`tvar'"
    ereturn local gvar "`gvar'"
    ereturn local outcome "`y'"
    ereturn local sourcecmd `"`sourcecmd'"'
    ereturn local jw_contract "scalar_joint_0.9.0"
    ereturn local jw_schema "scalar_v1"
    ereturn local jw_target "population_iw"
    ereturn local method "sunabraham_population"
    ereturn local jw_target_label "Population-share weighted mean over event times `postwindow'"
    ereturn local restriction_scale "cohort_outcome_contrast"
    ereturn local jw_restriction_scale "cohort_outcome_contrast"
    ereturn local cmd "jwselect_sa_iw"
    ereturn local cmdline `"jwselect_sa_iw `original0'"'
    ereturn local version "0.9.1"
    ereturn local predict "jwselect_p"
    ereturn local marginsnotok "_ALL"
    capture noisily signestimationsample `y' `ivar' `tvar' `gvar'
    local rc=_rc
    if `rc' quietly estimates restore `previous'
    quietly estimates drop `previous'
    if `rc' exit `rc'
    _jwselect_sa_iw_display
end

program define _jwselect_sa_iw_display
    version 16.1
    di as txt _n "Population interaction-weighted Sun-Abraham target"
    di as txt "Event-time window: " as res "`e(postwindow)'"
    di as txt "Estimate: " as res %10.5g _b[population_iw] ///
        as txt "    baseline SE: " as res %10.5g _se[population_iw]
    di as txt "Complete panel influence covariance includes estimated cohort shares."
    di as txt "Optional moments retain individual cohort/preperiod restrictions."
end
