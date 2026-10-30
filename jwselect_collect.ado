*! jwselect 0.7.0 23oct2026
program define jwselect_collect, rclass
    version 17.0
    syntax using/ [, PROGRESS]
    local directory `"`using'"'
    _jwselect_jobguard, directory(`"`directory'"') action(verify)
    _jwselect_jobguard, directory(`"`directory'"') action(unlocked)
    _jwselect_load
    tempname table diagnostics extrema output_map jobs settings
    local rngstate `c(rngstate)'
    capture noisily mata: jw_collect_stata()
    local rc = _rc
    quietly set rngstate `rngstate'
    if `rc' exit `rc'
    matrix colnames `table' = estimate se_fixed ll ul pvalue
    matrix colnames `diagnostics' = output objective nactive search_reps validation_reps ///
        evaluations starts_used converged_starts budget_exhausted search_value ///
        validated_value validation_gap radius_fraction
    matrix colnames `output_map' = kind cohort
    matrix colnames `jobs' = job_id output objective start
    matrix colnames `settings' = level gamma mcerr reps optstarts optmaxeval opttol optreps cachemb
    local k = colsof(`extrema')-2
    local hn
    if `k'>0 {
        forvalues j=1/`k' {
            local hn `hn' h`j'
        }
    }
    matrix colnames `extrema' = output objective `hn'
    return scalar N_jobs = `njobs'
    return matrix table = `table'
    return matrix diagnostics = `diagnostics'
    return matrix extrema = `extrema'
    return matrix output_map = `output_map'
    return matrix jobs = `jobs'
    return matrix settings = `settings'
    return local directory `"`directory'"'
    return local inference reselect
    di as result "Focused inference collected. Results are in r(table), r(diagnostics), and r(extrema)."
    di as text "Inspect search diagnostics; completion does not certify global extrema."
end
