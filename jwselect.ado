*! jwselect 0.4.0 23oct2026
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
        GAP(real 1) BASEline(string) ANTicipation(integer 0) ///
        SELECT(string) WBICScope(string) FOCUS(string) CANDidates(string) ///
        PREMatrix(name) CANDMatrix(name) JOINTMatrix(name) ///
        POSTWindow(numlist integer min=2 max=2) ///
        AGGregate(string) AGGWeights(name) VCE(string) ///
        INFERence(string) Level(cilevel) GAMma(real .01) ///
        MCErr(real .001) REPS(integer 19999) SEED(integer 12345) ///
        OPTStarts(integer 4) OPTMaxeval(integer 400) OPTTol(real .02) ///
        OPTReps(integer 1999) INFEROutputs(numlist integer) ///
        PREPare(string) CACHEMB(real 64) PROGRESS ///
        MAXCand(integer 4096) MAXJoint(integer 4096) NOCohort CELLEffects]

    if "`baseline'" == "" local baseline cohort
    if "`select'" == "" local select wbic
    local wbicscope_explicit = ("`wbicscope'" != "")
    if "`wbicscope'" == "" local wbicscope cohort
    if "`focus'" == "" local focus both
    if "`candidates'" == "" local candidates nested
    if "`aggregate'" == "" local aggregate equalcohort
    if "`inference'" == "" local inference auto
    foreach opt in baseline select wbicscope focus candidates aggregate vce inference {
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
    if !inlist("`inference'", "auto", "reselect", "envelope", "none") error 198
    if inlist("`inference'", "reselect", "envelope") & "`select'" != "focused" {
        di as err "inference(reselect or envelope) requires select(focused)"
        exit 198
    }
    if `cachemb'<0 | missing(`cachemb') error 198
    if "`inferoutputs'"!="" & "`select'"!="focused" {
        di as err "inferoutputs() requires select(focused)"
        exit 198
    }
    if "`inferoutputs'"!="" & !inlist("`inference'","none","auto","reselect") {
        di as err "inferoutputs() supports inference(none, auto, or reselect)"
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
    tempname rawb rawV rawcells rawtable
    tempvar sample
    local N_obs .
    local N_covariates .
    local N_clusters .
    local upstream_anticipation .
    local sourcecmd matrix
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
        if "`e(cmd)'" != "jwdid" {
            di as err "run jwdid first, or supply bmat(), vmat(), cells(), and nunits()"
            exit 301
        }
        quietly generate byte `sample' = e(sample)
        if "`vce'" == "" local vce unconditional
        _jwselect_jwdid, vce(`vce')
        matrix `rawb' = r(bgt)
        matrix `rawV' = r(Vgt)
        matrix `rawcells' = r(cells)
        matrix `rawtable' = r(upstream_table)
        local upstream_stat `r(upstream_stat)'
        local nunits = r(N_units)
        local N_obs = r(N_obs)
        local N_covariates = r(N_covariates)
        local N_clusters = r(N_clusters)
        local gap = r(gap)
        local upstream_anticipation = r(anticipation)
        local sourcecmd `"`r(sourcecmd)'"'
    }

    // Check the entry point as well as the completion/version marker.
    capture mata: assert(findexternal("jw_stata()") != NULL)
    local entryrc = _rc
    capture mata: assert(jw_version() == "0.4.0-joint1")
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
        capture mata: assert(jw_version() == "0.4.0-joint1")
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
        outz outzV outA outT outdiagnostics outextrema
    if "`select'" == "focused" & ("`inference'" != "none" | `"`prepare'"'!="") quietly set seed `seed'
    if "`select'" == "focused" & inlist("`inference'", "auto", "reselect") {
        di as txt "Computing focused reselection inference (numerical search and independent simulation)..."
    }
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
    if !`frommatrix' ereturn matrix upstream_table = `rawtable'
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
    ereturn scalar N_units = `nunits'
    ereturn scalar N_observations = `N_obs'
    ereturn scalar N_covariates = `N_covariates'
    ereturn scalar N_clusters = `N_clusters'
    ereturn scalar G = `_G'
    ereturn scalar K = `_K'
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
    ereturn local cmdline `"jwselect `0'"'
    ereturn local estat_cmd jwselect_estat
    ereturn local predict jwselect_p
    ereturn local version 0.4.0
    ereturn local build 0.4.0-joint1
    ereturn local cmd jwselect
    _jwselect_display, `nocohort' `celleffects'
end
