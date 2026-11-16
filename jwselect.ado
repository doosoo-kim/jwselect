*! jwselect 0.9.1 16nov2026
program define jwselect, eclass
    version 16.1
    if "`e(cmd)'" == "jwselect" {
        capture _jwselect_replay_options `0'
        if !_rc {
            _jwselect_display `0'
            exit
        }
    }
    tempname previous
    local held 0
    capture estimates store `previous'
    if !_rc local held 1
    local rngstate `c(rngstate)'
    capture noisily _jwselect_work `0'
    local rc = _rc
    quietly set rngstate `rngstate'
    if `rc' {
        if `held' quietly estimates restore `previous'
        else ereturn clear
    }
    if `held' quietly estimates drop `previous'
    if `rc' exit `rc'
end

program define _jwselect_replay_options
    version 16.1
    syntax [, NOCohort CELLEffects]
end

program define _jwselect_work, eclass
    version 16.1
    syntax [, BMAT(name) VMAT(name) CELLS(name) NUNits(real 0) ///
        GAP(real 1) BASEline(string) ANTicipation(numlist integer max=1 >=0) ///
        SELECT(string) WBICScope(string) FOCUS(string) CANDidates(string) ///
        PREMatrix(name) CANDMatrix(name) JOINTMatrix(name) ///
        POSTWindow(numlist integer min=2 max=2) ///
        AGGregate(string) AGGWeights(name) VCE(string) RESTRictions(string) ///
        INFERence(string) WBICMethod(string) MCMOde(string) MCFallback(string) ///
        WBICSolver(string) WBICEngine(string) MCBeta(numlist max=1 >0 <1) ///
        MCValidate(numlist integer max=1 >=100) MCSteps(numlist integer max=1 >=2) ///
        MCCells(numlist integer max=1 >=1) MCREFine MCRadius(numlist max=1 >=0) MCPMax(numlist max=1 >=0 <1) ///
        Level(cilevel) GAMma(real .01) MCErr(real .001) REPS(integer 19999) SEED(integer 12345) ///
        OPTStarts(integer 4) OPTMaxeval(integer 400) OPTTol(real .02) ///
        OPTReps(integer 1999) INFEROutputs(numlist integer) ///
        PREPare(string) CACHEMB(real 64) PROGRESS ///
        MAXCand(integer 4096) MAXJoint(integer 4096) NOCohort CELLEffects]

    local anticipation_explicit = ("`anticipation'" != "")
    if "`anticipation'" == "" local anticipation 0
    if "`baseline'" == "" local baseline cohort
    if "`select'" == "" local select wbic
    local wbicscope_explicit = ("`wbicscope'" != "")
    if "`wbicscope'" == "" local wbicscope cohort
    local focus_explicit = ("`focus'" != "")
    if "`focus'" == "" local focus both
    if "`candidates'" == "" local candidates nested
    local aggregate_explicit = ("`aggregate'" != "")
    if "`aggregate'" == "" local aggregate equalcohort
    if "`inference'" == "" local inference auto
    foreach opt in baseline select wbicscope focus candidates aggregate vce inference restrictions wbicmethod mcmode mcfallback wbicsolver wbicengine {
        local `opt' = lower("``opt''")
    }
    if !inlist("`baseline'", "cohort", "common") error 198
    if !inlist("`select'", "manual", "wbic", "focused") error 198
    if !inlist("`wbicscope'", "cohort", "joint") error 198
    if `wbicscope_explicit' & "`select'" != "wbic" {
        di as err "wbicscope() requires select(wbic)"
        exit 198
    }
    if !inlist("`focus'", "cohort", "aggregate", "both") error 198
    if !inlist("`candidates'", "nested", "all") error 198
    if !inlist("`aggregate'", "equalcohort", "equalcell") error 198
    if !inlist("`restrictions'", "", "index", "placebo") {
        di as err "restrictions() must be index (default) or placebo"
        exit 198
    }
    if !inlist("`inference'", "auto", "reselect", "envelope", "wbic", "none") error 198
    if inlist("`inference'", "reselect", "envelope") & "`select'" != "focused" {
        di as err "inference(reselect or envelope) requires select(focused)"
        exit 198
    }
    if "`inference'" == "wbic" & "`select'" != "wbic" {
        di as err "inference(wbic) requires select(wbic)"
        exit 198
    }
    if "`inference'"=="wbic" & !inlist("`wbicmethod'", "bootstrap", "mc") {
        di as err "inference(wbic) requires wbicmethod(bootstrap) or wbicmethod(mc)"
        exit 198
    }
    if "`wbicmethod'"!="" & "`inference'"!="wbic" {
        di as err "wbicmethod() requires inference(wbic)"
        exit 198
    }
    local bootstrap = ("`wbicmethod'"=="bootstrap")
    local mc = ("`wbicmethod'"=="mc")
    if "`wbicengine'"!="" & (!`bootstrap' | !inlist("`wbicengine'","integrated","direct")) {
        di as err "wbicengine(integrated or direct) requires wbicmethod(bootstrap)"
        exit 198
    }
    if `bootstrap' & "`wbicengine'"=="" local wbicengine integrated
    if "`mcmode'`wbicsolver'`mcfallback'`mcrefine'"!="" & !`mc' {
        di as err "McCloskey controls require inference(wbic) wbicmethod(mc)"
        exit 198
    }
    if !inlist("`mcmode'","","certified") | !inlist("`wbicsolver'","","grid") {
        di as err "McCloskey inference uses mcmode(certified) and wbicsolver(grid)"
        exit 198
    }
    if !inlist("`mcfallback'","","terminal","error") {
        di as err "mcfallback() must be terminal or error"
        exit 198
    }
    if `mc' & "`mcfallback'"=="" local mcfallback terminal
    local mccontrols mcbeta mcvalidate mcsteps mccells mcradius mcpmax
    local mcdefaults .4     9999       9       256     0        0
    local i 0
    foreach opt of local mccontrols {
        local ++i
        if "``opt''"!="" & !`mc' {
            di as err "`opt'() requires inference(wbic) wbicmethod(mc)"
            exit 198
        }
        if "``opt''"=="" local `opt' : word `i' of `mcdefaults'
        if missing(``opt'') error 198
    }
    if `mccells'>cond("`mcrefine'"=="",65536,4000000) {
        di as err "mccells() exceeds the budget for this McCloskey calculation"
        exit 198
    }
    if `cachemb'<0 | missing(`cachemb') error 198
    if "`inferoutputs'"!="" & "`select'"!="focused" & "`inference'"!="wbic" {
        di as err "inferoutputs() requires select(focused) or inference(wbic)"
        exit 198
    }
    if "`inferoutputs'"!="" & !inlist("`inference'","none","auto","reselect","wbic") {
        di as err "inferoutputs() supports inference(none, auto, reselect, or wbic)"
        exit 198
    }
    if `"`prepare'"'!="" {
        if "`select'"!="focused" | "`inference'"!="none" {
            di as err "prepare() requires select(focused) inference(none)"
            exit 198
        }
        capture mkdir `"`prepare'"'
        if _rc {
            di as err "prepare() requires a NEW directory with an existing parent"
            exit 602
        }
    }
    if `anticipation' < 0 | `gap' <= 0 | missing(`gap') error 198
    if `maxcand' < 1 | `maxjoint' < 1 error 198
    if `reps' < 1 | `gamma' <= 0 | `gamma' >= 1 error 198
    if `mcerr' <= 0 | `mcerr' >= 1 error 198
    if `seed' < 0 error 198
    if `optstarts' < 1 | `optmaxeval' < `optstarts' | `optreps' < 100 error 198
    if `opttol' <= 0 | `opttol' >= .5 | missing(`opttol') error 198
    if "`select'" != "manual" & "`prematrix'" != "" {
        di as err "prematrix() is available only with select(manual)"
        exit 198
    }
    if "`select'" == "manual" & "`candmatrix'`jointmatrix'" != "" {
        di as err "manual selection uses prematrix(), not candidate matrices"
        exit 198
    }
    if "`jointmatrix'" != "" & !( ("`select'" == "focused" & "`focus'" != "cohort") | ("`select'" == "wbic" & "`wbicscope'" == "joint") ) {
        di as err "jointmatrix() requires joint Wald-BIC or aggregate-focused selection"
        exit 198
    }
    foreach mat in bmat vmat cells prematrix candmatrix jointmatrix aggweights {
        if "``mat''" != "" confirm matrix ``mat''
    }
    local frommatrix = ("`bmat'`vmat'`cells'" != "")
    if `frommatrix' & ("`bmat'" == "" | "`vmat'" == "" | "`cells'" == "") {
        di as err "supply bmat(), vmat(), and cells() together"
        exit 198
    }
    if "`restrictions'" != "" {
        if `frommatrix' | "`e(cmd)'" != "jwdid" | !inlist("`e(cmd2)'", "logit", "poisson") {
            di as err "restrictions() requires a supported nonlinear jwdid logit or poisson fit"
            exit 198
        }
    }
    if !`frommatrix' & "`e(cmd)'" == "jwdid" & inlist("`e(cmd2)'", "logit", "poisson") {
        if "`restrictions'" == "" local restrictions index
    }
    local scalarinput = (!`frommatrix' & "`e(jw_schema)'" == "scalar_v1")
    if `scalarinput' {
        if "`postwindow'`aggweights'" != "" | "`baseline'" != "cohort" | `anticipation_explicit' | `aggregate_explicit' {
            di as err "Set the target, baseline, and population in the scalar first-stage command; they cannot be changed after fitting"
            exit 198
        }
        if `focus_explicit' & "`focus'" != "aggregate" {
            di as err "This source supplies one population target; use focus(aggregate)"
            exit 198
        }
        local focus aggregate
    }
    tempname rawd rawjointV rawbases rawtargetmap
    tempname rawb rawV rawcells rawtable
    local rawh
    local rawVh
    local rawVdh
    local rawhmap
    local index_terms
    local index_covariates
    local N_index 0
    tempvar sample
    local N_obs .
    local N_covariates .
    local N_clusters .
    local upstream_anticipation .
    local sourcecmd matrix
    local sourceestimator matrix
    local sourceversion
    local sourcemethod supplied
    local sourcemodel linear
    local restriction_scale outcome_contrast
    if `frommatrix' {
        if `nunits' <= 1 | missing(`nunits') {
            di as err "nunits() must give the independent-unit count, greater than one"
            exit 198
        }
        if "`vce'" == "" local vce supplied
        if "`vce'" != "supplied" {
            di as err "matrix input uses vce(supplied); its covariance convention is supplied by you"
            exit 198
        }
        matrix `rawb' = `bmat'
        matrix `rawV' = `vmat'
        matrix `rawcells' = `cells'
    }
    else {
        if !inlist("`e(cmd)'", "jwdid", "csdid", "jwselect_sa", "jwselect_rcs", "jwselect_sa_iw") {
            di as err "run a supported jwdid, csdid, or jwselect_sa fit first, or supply bmat(), vmat(), cells(), and nunits()"
            if "`e(cmd)'" == "reghdfe" di as err "For Sun-Abraham, fit with jwselect_sa so the complete event map and estimation sample can be verified"
            exit 301
        }
        local sourceestimator `e(cmd)'
        quietly generate byte `sample' = e(sample)
        if `scalarinput' {
            if "`vce'" == "" local vce supplied
            _jwselect_scalar, vce(`vce')
        }
        else if "`sourceestimator'" == "jwselect_sa" {
            if "`vce'" == "" local vce supplied
            _jwselect_sa, vce(`vce')
        }
        else if "`sourceestimator'" == "csdid" {
            if "`vce'" == "" local vce supplied
            _jwselect_csdid, vce(`vce')
        }
        else {
            if "`vce'" == "" local vce unconditional
            if "`restrictions'" == "index" {
                _jwselect_index, vce(`vce')
            }
            else _jwselect_jwdid, vce(`vce')
        }
        if `scalarinput' {
            matrix `rawd' = r(d)
            matrix `rawjointV' = r(V)
            matrix `rawbases' = r(bases)
            matrix `rawtargetmap' = r(targetmap)
            local scalar_target_label `"`r(target_label)'"'
            if r(K)>0 {
                tempname rawh rawhmap
                matrix `rawh' = r(h)
                matrix `rawhmap' = r(hmap)
            }
        }
        else {
            matrix `rawb' = r(bgt)
            matrix `rawV' = r(Vgt)
            matrix `rawcells' = r(cells)
            matrix `rawtable' = r(upstream_table)
        }
        local upstream_stat `r(upstream_stat)'
        local nunits = r(N_units)
        local N_obs = r(N_obs)
        local N_covariates = r(N_covariates)
        local N_clusters = r(N_clusters)
        local gap = r(gap)
        local upstream_anticipation = r(anticipation)
        if `scalarinput' local anticipation = `upstream_anticipation'
        local sourcecmd `"`r(sourcecmd)'"'
        local sourceversion `"`r(sourceversion)'"'
        local sourcemethod `"`r(method)'"'
        local sourcemodel `"`r(model)'"'
        local restriction_scale `"`r(restriction_scale)'"'
        if "`restriction_scale'" == "conditional_index" {
            local N_index = r(N_index)
            local upstream_vce `"`r(upstream_vce)'"'
            local index_covariance_factor = r(covariance_factor)
            local index_covariance_error = r(covariance_error)
            local index_margins_error = r(margins_error)
            local index_covariates `"`r(index_covariates)'"'
            if `N_index' > 0 {
                tempname rawh rawVh rawVdh rawhmap index_terms
                matrix `rawh' = r(h)
                matrix `rawVh' = r(Vh)
                matrix `rawVdh' = r(Vdh)
                matrix `rawhmap' = r(hmap)
                matrix `index_terms' = r(index_terms)
            }
        }
        if "`sourceestimator'" == "jwselect_sa" local sourceestimator eventstudyinteract
        // These backends supply baseline ATTs and zero-placebo moments.
        // Their effects cannot be rebased by subtracting another ATT.
        if !`scalarinput' & ("`sourceestimator'" == "csdid" | "`sourcemodel'" == "nonlinear") {
            if "`baseline'" != "cohort" {
                di as err "CS-DR and nonlinear adapters require baseline(cohort); alternate-baseline refits are not implemented"
                exit 198
            }
            local source_delta = `upstream_anticipation'
            if "`sourceestimator'" == "jwdid" local source_delta = `upstream_anticipation'-1
            if `anticipation_explicit' & `anticipation' != `source_delta' {
                di as err "anticipation() must match the source estimator's maintained baseline (`source_delta')"
                exit 198
            }
            local anticipation = `source_delta'
            tempname normbad
            forvalues j = 1/`=rowsof(`rawcells')' {
                if abs(`rawcells'[`j',2]-(`rawcells'[`j',1]-`gap'*(`source_delta'+1))) < 1e-8*`gap' {
                    mata: st_numscalar("`normbad'", max(abs(st_matrix("`rawV'")[`j',.])))
                    if `rawb'[1,`j'] != 0 | scalar(`normbad') != 0 {
                        di as err "The source baseline must be a structural zero with zero covariance"
                        exit 498
                    }
                }
            }
        }
    }

    // Check the entry point as well as the completion/version marker.
    capture mata: assert(findexternal("jw_stata()") != NULL)
    local entryrc = _rc
    capture mata: assert(jw_version() == "0.9.1" & jw_build() == "0.9.1-20261116")
    local versionrc = _rc
    if `entryrc' | `versionrc' {
        capture mata: assert(findexternal("jw_version()") != NULL)
        local markerrc = _rc
        if !`entryrc' | !`markerrc' {
            di as err "jwselect Mata code is incomplete or from a different build"
            di as err "Save your work, restart Stata, and rerun with the updated package"
            exit 499
        }
        quietly findfile _jwselect_mata.do
        local loader `"`r(fn)'"'
        // noisily is intentional: expose the FIRST compiler error.
        capture noisily do `"`loader'"'
        local loadrc = _rc
        if `loadrc' {
            di as err "jwselect Mata compilation failed; the first error is shown above"
            di as err "Loader: `loader'"
            exit `loadrc'
        }
        capture mata: assert(findexternal("jw_stata()") != NULL)
        local entryrc = _rc
        capture mata: assert(jw_version() == "0.9.1" & jw_build() == "0.9.1-20261116")
        local versionrc = _rc
        if `entryrc' | `versionrc' {
            di as err "The Mata loader ended without defining the complete jwselect module"
            di as err "Loader: `loader'"
            di as err "Run this file directly in a fresh Stata session and report its first error"
            exit 499
        }
    }
    tempname outb outV outtable outmap outselection outbases outpost outh ///
        outbcells outVcells outcellsel outweights outscores outfamily outjoint ///
        outz outzV outA outT outdiagnostics outextrema outblockmap outmcdiag outmccal outbbdiag
    if "`select'" == "focused" & ("`inference'" != "none" | `"`prepare'"'!="") quietly set seed `seed'
    if "`inference'" == "wbic" quietly set seed `seed'
    if "`select'" == "focused" & inlist("`inference'", "auto", "reselect") {
        di as txt "Computing focused reselection inference (numerical search and independent simulation)..."
    }
    if `bootstrap' di as txt "Computing centered WBIC bootstrap interval (`wbicengine', `reps' draws)..."
    if `mc' di as txt "Computing McCloskey intervals with continuum calibration and numerical padding..."
    tempname outpwdiag outgrid outprofile outterminal
    mata: jw_stata()
    if `"`prepare'"'!="" {
        _jwselect_jobguard, directory(`"`prepare'"') action(finalize)
    }

    local nobsopt
    if !missing(`N_obs') local nobsopt obs(`N_obs')
    local sampleopt
    if !`frommatrix' local sampleopt esample(`sample')
    tempname Vfixed ci pv
    matrix `Vfixed' = `outV'
    matrix `ci' = `outtable'[1...,3..4]
    matrix `pv' = `outtable'[1...,5]
    if "`select'" == "focused" ereturn post `outb', `sampleopt' `nobsopt'
    else ereturn post `outb' `outV', `sampleopt' `nobsopt'
    if !`frommatrix' & !`scalarinput' ereturn matrix upstream_table = `rawtable'
    ereturn matrix V_fixed = `Vfixed'
    ereturn matrix table = `outtable'
    ereturn matrix ci = `ci'
    ereturn matrix pvalue = `pv'
    ereturn matrix output_map = `outmap'
    ereturn matrix baselines = `outbases'
    ereturn matrix postmap = `outpost'
    if `_K' > 0 {
        ereturn matrix hmap = `outh'
        ereturn matrix selection = `outselection'
        ereturn matrix cell_selection = `outcellsel'
        if ("`select'" == "focused" & "`focus'" != "cohort") | ("`select'" == "wbic" & "`wbicscope'" == "joint") ereturn matrix joint_candidates = `outjoint'
    }
    ereturn matrix b_cells = `outbcells'
    ereturn matrix V_cells = `outVcells'
    ereturn matrix weights = `outweights'
    ereturn matrix scores = `outscores'
    ereturn matrix candidates = `outfamily'
    ereturn matrix z = `outz'
    ereturn matrix V_z = `outzV'
    ereturn matrix loading = `outA'
    ereturn matrix targets = `outT'
    if "`select'" == "focused" & inlist("`inference'", "auto", "reselect") {
        ereturn matrix inference_diagnostics = `outdiagnostics'
        capture confirm matrix `outextrema'
        if !_rc ereturn matrix inference_extrema = `outextrema'
    }
    if `bootstrap' {
        ereturn local wbicmethod bootstrap
        ereturn matrix wbic_pointwise_diagnostics = `outpwdiag'
        ereturn local wbicengine `wbicengine'
        ereturn local inference_bound centered_reselection_t
    }
    if `mc' {
        ereturn local wbicmethod mc
        ereturn local wbicsolver grid
        ereturn local mcmode certified
        ereturn local mcfallback `mcfallback'
        ereturn local inference_bound mc_certified_calibration
        ereturn matrix mc_diagnostics = `outmcdiag'
        ereturn matrix mc_profile_diagnostics = `outprofile'
        ereturn matrix mc_terminal_diagnostics = `outterminal'
        ereturn local mc_unbounded_rule "unbounded=1 denotes [-inf,+inf]; endpoints stored missing"
        capture confirm matrix `outmccal'
        if !_rc ereturn matrix mc_calibration = `outmccal'
        capture confirm matrix `outbbdiag'
        if !_rc ereturn matrix bb_diagnostics = `outbbdiag'
        capture confirm matrix `outgrid'
        if !_rc ereturn matrix wbic_bb_grid_diagnostics = `outgrid'
        ereturn scalar mc_domain_radius = `_mc_domain_radius'
        ereturn scalar mc_exterior_radius = `_mc_exterior_radius'
        ereturn scalar mcrefine = ("`mcrefine'"!="")
        foreach control of local mccontrols {
            ereturn scalar `control' = ``control''
        }
    }
    ereturn scalar scalar_target = `scalarinput'
    if `scalarinput' ereturn local target_label `"`scalar_target_label'"'
    ereturn scalar N_units = `nunits' 
    ereturn scalar N_observations = `N_obs'
    ereturn scalar N_covariates = `N_covariates'
    ereturn scalar N_clusters = `N_clusters'
    ereturn scalar G = `_G'
    ereturn scalar K = `_K'
    ereturn scalar K_blocks = `_K_blocks'
    if "`restriction_scale'" == "conditional_index" {
        ereturn local upstream_vce `upstream_vce'
        ereturn local index_basis raw_covariates
        ereturn scalar index_covariance_factor = `index_covariance_factor'
        ereturn scalar index_covariance_error = `index_covariance_error'
        ereturn scalar index_margins_error = `index_margins_error'
    }
    if "`restriction_scale'" == "conditional_index" & `_K' > 0 {
        ereturn matrix blockmap = `outblockmap'
        ereturn matrix index_b = `rawh'
        ereturn matrix V_index = `rawVh'
        ereturn matrix Cov_att_index = `rawVdh'
        ereturn matrix index_hmap = `rawhmap'
        ereturn matrix index_terms = `index_terms'
        ereturn matrix upstream_cells = `rawcells'
        ereturn local index_covariates `"`index_covariates'"'
        ereturn local index_basis raw_covariates
    }
    ereturn scalar P = `_P'
    ereturn scalar gap = `gap'
    ereturn scalar anticipation = `anticipation'
    ereturn scalar upstream_anticipation = `upstream_anticipation'
    ereturn scalar level = `level'
    ereturn scalar gamma = `gamma'
    ereturn scalar mcerr = `mcerr'
    ereturn scalar reps = `reps'
    ereturn scalar seed = `seed'
    ereturn scalar maxcand = `maxcand'
    ereturn scalar maxjoint = `maxjoint'
    ereturn scalar cachemb = `cachemb'
    ereturn local inferoutputs `inferoutputs'
    ereturn local prepared_directory `"`prepare'"'
    ereturn scalar optstarts = `optstarts'
    ereturn scalar optmaxeval = `optmaxeval'
    ereturn scalar opttol = `opttol'
    ereturn scalar optreps = `optreps'
    ereturn local upstream_stat `upstream_stat'
    ereturn local baseline `baseline'
    ereturn local selection_method `select'
    if "`select'" == "wbic" ereturn local wbic_scope `wbicscope'
    ereturn local focus `focus'
    ereturn local candidate_type `candidates'
    if "`candmatrix'" != "" ereturn local candidate_type custom
    ereturn local aggregate `aggregate'
    if "`aggweights'" != "" ereturn local aggregate custom
    ereturn local covariance `vce'
    ereturn local inference `inference'
    if "`inference'" == "auto" {
        if "`select'" == "focused" ereturn local inference reselect
        if "`select'" == "wbic" ereturn local inference pointwise_oracle
        if "`select'" == "manual" ereturn local inference fixed_restrictions
    }
    ereturn local postwindow `postwindow'
    ereturn local sourcecmd `"`sourcecmd'"'
    ereturn local source_estimator `sourceestimator'
    ereturn local source_version `sourceversion'
    ereturn local source_method `sourcemethod'
    ereturn local source_model `sourcemodel'
    ereturn local restriction_scale `restriction_scale'
    if "`sourcemodel'" == "nonlinear" {
        if "`restrictions'" == "" local restrictions index
        ereturn local restrictions `restrictions'
    }
    ereturn local cmdline `"jwselect `0'"'
    ereturn local estat_cmd jwselect_estat
    ereturn local predict jwselect_p
    ereturn local version 0.9.1
    ereturn local build 0.9.1-20261116
    ereturn local cmd jwselect
    _jwselect_display, `nocohort' `celleffects'
end
