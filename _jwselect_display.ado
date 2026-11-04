*! jwselect display 0.9.0 04nov2026
program define _jwselect_display
    version 16.1
    syntax [, NOCohort CELLEffects]
    if "`e(cmd)'" != "jwselect" error 301
    di as txt _n "jwselect: preperiod selection and covariance-weighted moment adjustment"
    local scalarinput = (e(scalar_target)==1)
    if `scalarinput' di as txt "Population target: " as res `"`e(target_label)'"'
    local source `e(source_estimator)' 
    if "`source'" == "" local source jwdid
    local base "Cohort-specific (g-1)"
    if "`e(baseline)'" == "common" local base "Common (q-1)"
    if e(anticipation)>0 {
        local shift = 1+e(anticipation)
        local base "Cohort-specific (g-`shift')"
        if "`e(baseline)'" == "common" local base "Common (q-`shift')"
    }
    if e(gap)!=1 {
        local base "Cohort-specific: g-gap*(1+delta)"
        if "`e(baseline)'" == "common" local base "Common: q-gap*(1+delta)"
    }
    if `scalarinput' local base "Fixed by source; see baselines below"
    local selection "Cohort-wise Wald-BIC"
    if "`e(selection_method)'" == "wbic" & "`e(wbic_scope)'" == "joint" local selection "Joint Wald-BIC across cohorts"
    if "`e(selection_method)'" == "manual" local selection "Manually specified preperiods"
    if "`e(selection_method)'" == "focused" local selection "Focused selection"
    local family "Nested sets: recent to earlier"
    if "`e(candidate_type)'" == "all" local family "All subsets of eligible preperiods"
    if "`e(candidate_type)'" == "custom" local family "User-specified candidate sets"
    if "`e(wbic_scope)'" == "joint" local family "Joint configurations; see estat selection"
    if "`e(selection_method)'" == "manual" local family "Not searched (manual selection)"
    local agg "Uniform weights across cohorts"
    if "`e(aggregate)'" == "equalcell" local agg "Uniform weights across post cells"
    if "`e(aggregate)'" == "custom" local agg "User-specified fixed cell weights"
    if `scalarinput' local agg "Population target fixed in source fit"
    local inf "Pointwise oracle normal inference"
    if "`e(inference)'" == "fixed_restrictions" local inf "Normal; fixed restrictions"
    if "`e(inference)'" == "envelope" local inf "Conservative Gaussian envelope"
    if "`e(inference)'" == "reselect" local inf "Two-step numerical reselection"
    if "`e(inference)'" == "wbic" local inf "Selection-adjusted WBIC (`e(wbicmethod)' bound)"
    if "`e(wbicmethod)'" == "mc" local inf "McCloskey size calibration"
    if "`e(wbicmethod)'" == "bootstrap" local inf "Centered WBIC bootstrap-t (`e(wbicengine)')"
    if "`e(inference)'" == "none" local inf "CIs and p-values disabled"
    di as txt _n "Research design and method"
    di as txt "{hline 38}{c TT}{hline 39}"
    di as txt "Source estimator" _col(39) "{c |} " as res "`source' `e(source_method)'"
    di as txt "Independent units" _col(39) "{c |} " as res %9.0g e(N_units)
    if !missing(e(N_observations)) di as txt "Observations" _col(39) "{c |} " as res %9.0g e(N_observations)
    if !missing(e(N_clusters)) di as txt "Clusters" _col(39) "{c |} " as res %9.0g e(N_clusters)
    di as txt "Treatment cohorts" _col(39) "{c |} " as res %9.0g e(G)
    if !missing(e(N_covariates)) di as txt "Covariates" _col(39) "{c |} " as res %9.0g e(N_covariates)
    if `scalarinput' di as txt "Population targets" _col(39) "{c |} " as res %9.0g e(P)
    else di as txt "Target post-treatment cells" _col(39) "{c |} " as res %9.0g e(P)
    if "`e(restriction_scale)'" == "conditional_index" {
        di as txt "Eligible preperiod blocks" _col(39) "{c |} " as res %9.0g e(K_blocks)
        di as txt "Index-coefficient restrictions" _col(39) "{c |} " as res %9.0g e(K)
    }
    else di as txt "Eligible over-identifying preperiods" _col(39) "{c |} " as res %9.0g e(K)
    di as txt "{hline 38}{c +}{hline 39}"
    di as txt "Baseline" _col(39) "{c |} " as res "`base'"
    di as txt "Allowed anticipation (delta)" _col(39) "{c |} " as res %9.0g e(anticipation)
    di as txt "Period spacing (gap)" _col(39) "{c |} " as res %9.0g e(gap)
    di as txt "Selection method" _col(39) "{c |} " as res "`selection'"
    di as txt "Candidate preperiod sets" _col(39) "{c |} " as res "`family'"
    di as txt "Aggregation" _col(39) "{c |} " as res "`agg'"
    if "`e(selection_method)'" == "focused" {
        local focus "Cohort-wise and aggregate-focused"
        if "`e(focus)'" == "cohort" local focus "Cohort-wise focused"
        if "`e(focus)'" == "aggregate" local focus "Aggregate-focused joint selection"
        di as txt "Selection target" _col(39) "{c |} " as res "`focus'"
    }
    di as txt "Covariance" _col(39) "{c |} " as res "`e(covariance)'"
    di as txt "Inference" _col(39) "{c |} " as res "`inf'"
    local window "All available post-treatment periods"
    if "`e(postwindow)'" != "" local window "Event times `e(postwindow)' (inclusive)"
    if `scalarinput' local window "Fixed by source target definition"
    di as txt "Post-treatment window" _col(39) "{c |} " as res "`window'"
    di as txt "{hline 38}{c BT}{hline 39}"
    di as txt "Preperiod counts are cohort-period pairs, excluding maintained baselines."
    if !`scalarinput' di as txt "Each cohort effect averages its included post-treatment cells equally."
    else di as txt "All candidate restrictions use the same prespecified population target and its full covariance."
    if "`e(restriction_scale)'" == "outcome_placebo" {
        di as txt "Optional restrictions set average outcome-scale placebo effects to zero."
        di as txt "J equals Wald for this moment-adjustment criterion; no upstream model is refitted."
    }
    if "`e(restriction_scale)'" == "conditional_index" {
        di as txt "Each optional date selects its whole conditional-index coefficient block."
        di as txt "Wald-BIC counts independent coefficient restrictions, not selected dates."
        di as txt "ATT estimates use target-level MD with the full joint covariance."
        di as txt "This is not a constrained nonlinear refit; adjusted ATTs need not preserve outcome bounds."
    }
    if "`e(baseline)'" == "common" di as txt "q denotes the earliest treatment cohort."
    tempname upstream
    capture matrix `upstream' = e(upstream_table)
    local matrices : e(matrices)
    if strpos(" `matrices' ", " upstream_table ") {
        di as txt _n "Initial cohort-time estimates (from `source')"
        di as txt "Time is calendar time. Zero-SE reference rows are normalizations."
        if "`source'" == "csdid" di as txt "Pointwise normal inference below uses csdid's analytical covariance, before selection."
        else di as txt "Inference below is copied from `source'; it precedes preperiod selection."
        if "`e(upstream_vce)'" == "conditional" {
            di as txt "This source comparison table uses conditional margins SEs."
            di as txt "The adjusted results below use vce(`e(covariance)') and the full joint covariance."
        }
        local statlabel "`e(upstream_stat)'"
        if !inlist("`statlabel'", "t", "z") local statlabel "b/se"
        di as txt %7s "Cohort" " " %7s "Time" " {c |} " ///
            %10s "Estimate" " " %10s "SE" " " %8s "`statlabel'" " " ///
            %10s "p-value" " " %10s "CI lower" " " %10s "CI upper"
        di as txt "{hline 16}{c +}{hline 66}"
        forvalues i=1/`=rowsof(`upstream')' {
            di as res %7.0g `upstream'[`i',1] " " %7.0g `upstream'[`i',2] ///
                as txt " {c |} " as res %10.5g `upstream'[`i',3] " " ///
                %10.5g `upstream'[`i',4] " " %8.3f `upstream'[`i',5] " " ///
                %10.4g `upstream'[`i',6] " " %10.5g `upstream'[`i',7] " " ///
                %10.5g `upstream'[`i',8]
        }
        di as txt "{hline 16}{c BT}{hline 66}"
    }
    di as txt _n "Treatment effects after preperiod selection"
    if "`e(selection_method)'" == "focused" {
        di as txt "SEs are fixed-selection diagnostics."
        if "`e(inference)'" == "none" di as txt "Selection-adjusted inference was disabled."
        else if "`e(inference)'" == "envelope" di as txt "CIs and p-values use the conservative Gaussian envelope."
        else {
            di as txt "CIs and p-values repeat focused selection in Gaussian draws."
            di as txt "The numerical nuisance search is not a certified global bound; see estat inference."
            tempname diagnostics
            matrix `diagnostics' = e(inference_diagnostics)
            local exhausted 0
            forvalues i=1/`=rowsof(`diagnostics')' {
                if `diagnostics'[`i',9] == 1 local ++exhausted
            }
            if `exhausted' > 0 di as txt "Search budget exhausted for `exhausted' objectives; inspect estat inference."
        }
    }
    else if "`e(inference)'" == "wbic" {
        if "`e(wbicmethod)'"=="mc" {
            tempname mcterminals
            mata: st_numscalar("`mcterminals'",sum(st_matrix("e(mc_diagnostics)")[,14]:==7))
            if scalar(`mcterminals')>0 {
                di as txt "Analytic McCloskey terminal used for " scalar(`mcterminals') " output(s)."
                di as txt "Inspect the terminal normal score and width with estat inference."
            }
            di as txt "McCloskey calibration reports CIs; calibrated p-values are not supplied."
        }
        else if "`e(wbicmethod)'"=="bootstrap" di as txt "Studentized reselection about the observed WBIC estimate. Coverage justification is pointwise."
        di as txt "Reported SEs and e(V) are fixed-selection diagnostics. Use e(ci) for interval inference."
        di as txt "Intervals are marginal. Gaussian approximation and regularity conditions apply."
        di as txt "See estat inference for method-specific diagnostics."
    }
    else if "`e(selection_method)'" == "wbic" & "`e(inference)'" != "none" ///
        di as txt "Normal inference uses the pointwise oracle approximation after Wald-BIC."
    else if "`e(inference)'" == "none" di as txt "Confidence intervals and p-values were disabled."
    tempname shown
    matrix `shown' = e(table)
    if "`nocohort'" != "" & !`scalarinput' matrix `shown' = `shown'[`=e(G)+1'...,1...]
    matlist `shown', format(%11.5g) rowtitle(Target)
    if "`e(wbicmethod)'"=="mc" {
        forvalues i=1/`=rowsof(e(mc_diagnostics))' {
            if el(e(mc_diagnostics),`i',14)==8 {
                if "`nocohort'"=="" | `i'>e(G) ///
                    di as txt "Target row `i': CI [-inf,+inf] (whole real line; arithmetic or regularity exception)."
            }
        }
    }
    di as txt _n "Cohort baselines and selected over-identifying preperiods:"
    _jwselect_periods
    if "`e(selection_method)'" == "focused" & "`e(focus)'" == "both" {
        di as txt _n "Aggregate-focused over-identifying preperiods:"
        _jwselect_periods, aggregate
    }
    if "`e(aggregate)'" == "custom" di as txt "Custom weights apply to cells; the aggregate need not average displayed cohort summaries."
    if "`celleffects'" != "" jwselect_estat cells
end

program define _jwselect_periods
    version 16.1
    syntax [, AGGregate]
    tempname base hm sel
    matrix `base' = e(baselines)
    if e(K) > 0 {
        matrix `hm' = e(hmap)
        matrix `sel' = e(selection)
    }
    forvalues g=1/`=e(G)' {
        local cohort = `base'[`g',1]
        local baseperiod = `base'[`g',2]
        local row `g'
        if e(scalar_target)==1 local row 1
        if "`aggregate'" != "" local row = colsof(e(b))
        local periods
        if e(K) > 0 {
            forvalues h=1/`=e(K)' {
                if `hm'[`h',1] == `cohort' & `sel'[`row',`h'] == 1 {
                    local date = `hm'[`h',2]
                    local periods `periods' `date'
                }
            }
        }
        local periods : list uniq periods
        if "`periods'" == "" local periods "(none)"
        di as txt "  cohort " as res %10.0g `cohort' ///
            as txt "  baseline " as res %10.0g `baseperiod' ///
            as txt "  over-identifying: " as res "`periods'"
    }
    if "`e(selection_method)'" == "wbic" & "`e(wbic_scope)'" == "joint" ///
        di as txt "  The final joint MD fit uses the joint Wald-BIC selected restrictions."
    else if inlist("`e(selection_method)'", "manual", "wbic") ///
        di as txt "  The final joint MD fit uses the union of these cohort-specific restrictions."
    if "`e(selection_method)'" == "focused" & "`e(focus)'" == "aggregate" {
        if e(scalar_target)==1 di as txt "  The population estimate uses the jointly selected target-focused restrictions."
        else di as txt "  Cohort estimates use the jointly selected aggregate-focused restrictions."
    }
end
