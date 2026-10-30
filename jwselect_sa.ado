*! jwselect 0.7.0 30oct2026: validated Sun-Abraham fitting wrapper
program define jwselect_sa, eclass
    version 16.1
    if replay() {
        if `"`0'"' != "" {
            di as err "jwselect_sa replay accepts no options; restore the source and run jwselect for selection"
            exit 198
        }
        if "`e(cmd)'" != "jwselect_sa" error 301
        _jwselect_sa_display
        exit
    }
    local original0 `"`0'"'
    syntax varname(numeric) [if] [in], IVAR(varname numeric) ///
        TVAR(varname numeric) GVAR(varname numeric) [Level(cilevel)]
    marksample touse
    markout `touse' `ivar' `tvar' `gvar'
    tempname previous
    local held 0
    capture estimates store `previous'
    if !_rc local held 1
    local rngstate `c(rngstate)'
    local morestate `c(more)'
    preserve
    capture noisily _jwselect_sa_fit `varlist', touse(`touse') ///
        ivar(`ivar') tvar(`tvar') gvar(`gvar') level(`level')
    local rc = _rc
    restore
    quietly set rngstate `rngstate'
    quietly set more `morestate'
    if !`rc' {
        // Rebuild e(sample) in the original data order, after restore. Signing
        // before restoring order would invalidate the signature immediately.
        capture noisily ereturn repost, esample(`touse')
        local rc = _rc
        if !`rc' {
            ereturn local cmdline `"jwselect_sa `original0'"'
            capture noisily signestimationsample `varlist' `ivar' `tvar' `gvar'
            local rc = _rc
        }
    }
    if `rc' {
        if `held' quietly estimates restore `previous'
        else ereturn clear
    }
    if `held' quietly estimates drop `previous'
    if `rc' exit `rc'
    _jwselect_sa_display
end

program define _jwselect_sa_fit, eclass
    version 16.1
    syntax varname(numeric), TOUSE(varname) IVAR(varname numeric) ///
        TVAR(varname numeric) GVAR(varname numeric) Level(real)
    local y `varlist'
    if "`y'" == "`ivar'" | "`y'" == "`tvar'" | "`y'" == "`gvar'" | ///
        "`ivar'" == "`tvar'" | "`ivar'" == "`gvar'" | "`tvar'" == "`gvar'" {
        di as err "Outcome, panel identifier, calendar time, and treatment cohort must be different variables"
        exit 198
    }
    foreach command in eventstudyinteract reghdfe avar ftools {
        capture which `command'
        if _rc {
            di as err "jwselect_sa requires the optional `command' package; see help jwselect_sa for installation"
            exit 499
        }
    }
    // This adapter checks the public arrays and the underlying coefficient
    // layout of the author's 0.1 release. Unknown layouts fail explicitly.
    quietly findfile eventstudyinteract.ado
    local sourcefile `"`r(fn)'"'
    tempname fh
    file open `fh' using `"`sourcefile'"', read text
    file read `fh' sourceheader
    file close `fh'
    if !regexm(lower(`"`sourceheader'"'), "version[ ]+0[.]1[ ]+24jan2022") {
        di as err "jwselect_sa is validated against eventstudyinteract 0.1 (24jan2022)"
        di as err "The installed source has an unrecognized version; use the documented source or a validated matrix interface"
        exit 198
    }
    quietly count if `touse'
    local nobs = r(N)
    if `nobs' == 0 error 2000
    capture assert `tvar'>0 & `tvar'==floor(`tvar') & ///
        `gvar'>=0 & `gvar'==floor(`gvar') if `touse'
    if _rc {
        di as err "jwselect_sa requires positive integer calendar times and nonnegative integer treatment cohorts; never treated is 0"
        exit 198
    }
    tempvar ncell ni gmin gmax idtag never cohort event
    quietly bysort `ivar' `tvar': egen long `ncell' = total(`touse')
    capture assert `ncell'==1 if `touse'
    if _rc {
        di as err "The estimation sample must contain one observation per panel unit and calendar period"
        exit 459
    }
    quietly bysort `ivar': egen long `ni' = total(`touse')
    quietly by `ivar': egen double `gmin' = min(cond(`touse',`gvar',.))
    quietly by `ivar': egen double `gmax' = max(cond(`touse',`gvar',.))
    capture assert `gmin'==`gmax' if `touse'
    if _rc {
        di as err "Treatment cohort must be constant within each sampled panel unit"
        exit 459
    }
    quietly levelsof `tvar' if `touse', local(times)
    local nt : word count `times'
    local first : word 1 of `times'
    local last : word `nt' of `times'
    if `nt'<2 | `nt'!=`last'-`first'+1 {
        di as err "jwselect_sa requires at least two consecutive, unit-spaced calendar periods"
        exit 198
    }
    capture assert `ni'==`nt' if `touse'
    if _rc {
        di as err "jwselect_sa requires a balanced estimation sample with common calendar-time support"
        exit 459
    }
    quietly egen byte `idtag' = tag(`ivar') if `touse'
    quietly count if `idtag'==1
    local nu = r(N)
    quietly count if `idtag'==1 & `gvar'==0
    if r(N)<2 {
        di as err "jwselect_sa requires at least two observed never-treated panel units coded gvar=0"
        exit 198
    }
    capture assert `gvar'==0 | (`gvar'>`first' & `gvar'<=`last') if `touse'
    if _rc {
        di as err "Every treated cohort must have an observed g-1 baseline and at least one observed treated period"
        di as err "Always-treated units and cohorts outside the observed time support are not supported"
        exit 198
    }
    quietly levelsof `gvar' if `touse' & `gvar'>0, local(groups)
    local ng : word count `groups'
    if `ng'<1 {
        di as err "At least one treated cohort is required"
        exit 198
    }
    foreach gg of local groups {
        quietly count if `idtag'==1 & `gvar'==`gg'
        if r(N)<2 {
            di as err "Each treatment cohort must contain at least two sampled panel units"
            exit 198
        }
    }
    // Set excluded observations to control: the author's cohort enumeration
    // is outside its touse condition. This keeps if/in samples self-contained.
    quietly generate byte `never' = (`gvar'==0 | !`touse')
    quietly generate double `cohort' = `gvar' if `touse' & `gvar'>0
    quietly generate double `event' = `tvar'-`gvar' if `touse' & `gvar'>0
    quietly levelsof `event' if `touse' & `gvar'>0 & `event'!=-1, local(events)
    local relvars
    foreach hh of local events {
        tempvar rel
        quietly generate byte `rel' = (`event'==`hh') if `touse'
        quietly replace `rel' = 0 if !`touse'
        local relvars `relvars' `rel'
    }
    local nr : word count `events'
    quietly eventstudyinteract `y' `relvars' if `touse', ///
        cohort(`cohort') control_cohort(`never') absorb(`ivar' `tvar') vce(cluster `ivar')
    if "`e(cmd)'"!="reghdfe" | "`e(vce)'"!="cluster" | ///
        e(N)!=`nobs' | e(N_clust)!=`nu' | missing(e(df_r)) | e(df_r)<=0 {
        di as err "eventstudyinteract did not return the expected unit-clustered full regression and sample"
        exit 498
    }
    capture assert e(sample)==`touse'
    if _rc {
        di as err "The Sun-Abraham regression changed the validated balanced estimation sample"
        exit 498
    }
    tempname B0 V0 BI DI B V C L VI EV miss bad factor pdratio
    foreach name in b V b_interact V_interact {
        capture confirm matrix e(`name')
        if _rc {
            di as err "eventstudyinteract did not return e(`name')"
            exit 498
        }
    }
    matrix `B0'=e(b)
    matrix `V0'=e(V)
    matrix `BI'=e(b_interact)
    matrix `DI'=e(V_interact)
    local k=colsof(`B0')
    if rowsof(`B0')!=1 | `k'!=`ng'*`nr'+1 | ///
        rowsof(`V0')!=`k' | colsof(`V0')!=`k' | ///
        rowsof(`BI')!=`ng' | colsof(`BI')!=`nr' | ///
        rowsof(`DI')!=`ng' | colsof(`DI')!=`nr' {
        di as err "Unrecognized eventstudyinteract coefficient/covariance layout"
        exit 498
    }
    local brows : rownames `BI'
    local drows : rownames `DI'
    local bcols : colnames `BI'
    local dcols : colnames `DI'
    local fullnames : colnames `B0'
    local vcols : colnames `V0'
    local vrows : rownames `V0'
    local tail : word `k' of `fullnames'
    if "`brows'"!="`groups'" | "`drows'"!="`groups'" | ///
        "`bcols'"!="`relvars'" | "`dcols'"!="`relvars'" | ///
        "`tail'"!="_cons" | "`fullnames'"!="`vcols'" | "`fullnames'"!="`vrows'" {
        di as err "Sun-Abraham cohort, event, or covariance labels do not match the fitting design"
        exit 498
    }
    mata: st_numscalar("`miss'", hasmissing(st_matrix("`B0'")) | hasmissing(st_matrix("`V0'")) | hasmissing(st_matrix("`BI'")) | hasmissing(st_matrix("`DI'")))
    if scalar(`miss') {
        di as err "Sun-Abraham returned missing coefficient or covariance entries"
        exit 498
    }
    local m=`ng'*`nt'
    matrix `C'=J(`m',2,.)
    matrix `L'=J(`m',`k',0)
    local row 0
    local gi 0
    foreach gg of local groups {
        local ++gi
        foreach tt of local times {
            local ++row
            matrix `C'[`row',1]=`gg'
            matrix `C'[`row',2]=`tt'
            local hh=`tt'-`gg'
            if `hh'!=-1 {
                local ei : list posof "`hh'" in events
                if !`ei' {
                    di as err "A supported cohort-time cell is missing from the event design"
                    exit 498
                }
                local j=(`ei'-1)*`ng'+`gi'
                if abs(`B0'[1,`j']-`BI'[`gi',`ei'])>1e-10*max(1,abs(`BI'[`gi',`ei'])) | ///
                    abs(`V0'[`j',`j']-`DI'[`gi',`ei'])>1e-10*max(1,abs(`DI'[`gi',`ei'])) | ///
                    `V0'[`j',`j']<=0 {
                    di as err "A supported SA cell is omitted, unidentified, or inconsistent with the underlying regression"
                    exit 498
                }
                matrix `L'[`row',`j']=1
            }
        }
    }
    matrix `B'=`B0'*`L''
    matrix `V'=`L'*`V0'*`L''
    // Independent cell-mean and iid-panel influence check (v46_12, SA IF).
    // e(sa_Vif) has no regression degrees-of-freedom correction; the adapter
    // uses e(sa_Vgt), retaining the supplied reghdfe covariance convention.
    quietly sort `ivar' `tvar'
    tempname X Y G IF DD TG CG DTE DCE BE IX VF SC FA ER VV EE vifpos
    mata: `X'=st_data(.,("`y'","`gvar'"),"`touse'")
    mata: `Y'=colshape(`X'[,1],`nt')
    mata: `G'=`X'[(0::(`nu'-1)):*`nt':+1,2]
    mata: `IF'=J(`nu',`m',0)
    mata: `BE'=J(1,`m',0)
    mata: `CG'=selectindex(`G':==0)
    forvalues j=1/`m' {
        local gg=`C'[`j',1]
        local ix=`C'[`j',2]-`first'+1
        local bse=`gg'-`first'
        mata: `TG'=selectindex(`G':==`gg')
        mata: `DD'=`Y'[,`ix']-`Y'[,`bse']
        mata: `DTE'=mean(`DD'[`TG'])
        mata: `DCE'=mean(`DD'[`CG'])
        mata: `BE'[1,`j']=`DTE'-`DCE'
        mata: `IF'[`TG',`j']=(`nu'/length(`TG')):*(`DD'[`TG']:-`DTE')
        mata: `IF'[`CG',`j']=-(`nu'/length(`CG')):*(`DD'[`CG']:-`DCE')
    }
    mata: `VF'=quadcross(`IF',`IF')/`nu'^2
    mata: `SC'=max((1,max(abs(`Y'))))
    mata: `ER'=max(abs(`BE'-st_matrix("`B'")))/`SC'
    mata: `IX'=selectindex(diagonal(`VF'):>0)
    mata: `FA'=.
    mata: `VV'=st_matrix("`V'")
    mata: st_numscalar("`vifpos'",length(`IX'))
    if scalar(`vifpos')>0 {
        mata: `FA'=`VV'[`IX'[1],`IX'[1]]/`VF'[`IX'[1],`IX'[1]]
    }
    mata: `ER'=max((`ER',max(abs(st_matrix("`V'")-`FA':*`VF'))/max((1e-12,max(abs(st_matrix("`V'")))))))
    mata: `IX'=selectindex(st_matrix("`C'")[,2]:!=st_matrix("`C'")[,1]:-1)
    mata: `EE'=symeigenvalues(`VV'[`IX',`IX'])
    mata: st_numscalar("`pdratio'",min(`EE')/max(abs(`EE')))
    mata: st_numscalar("`bad'",`ER')
    mata: st_numscalar("`factor'",`FA')
    mata: st_matrix("`VI'",`VF')
    mata: mata drop `X' `Y' `G' `IF' `DD' `TG' `CG' `DTE' `DCE' `BE' `IX' `VF' `SC' `FA' `ER' `VV' `EE'
    if missing(scalar(`bad'),scalar(`factor')) | scalar(`bad')>1e-6 | scalar(`factor')<=0 {
        di as err "The SA regression does not match the complete cell-mean contrasts and joint panel influence covariance"
        exit 498
    }
    if missing(scalar(`pdratio')) | scalar(`pdratio')<=1e-10 {
        di as err "The identified Sun-Abraham cells require a positive-definite joint covariance"
        di as err "The fitted covariance is singular or numerically ill-conditioned; more independent units or a separately justified design are required"
        exit 498
    }
    matrix colnames `C'=cohort time
    local names
    forvalues j=1/`m' {
        local names `names' sa`j'
    }
    matrix colnames `B'=`names'
    matrix rownames `V'=`names'
    matrix colnames `V'=`names'
    matrix rownames `VI'=`names'
    matrix colnames `VI'=`names'
    matrix `EV'=J(1,`nr',.)
    forvalues j=1/`nr' {
        local val : word `j' of `events'
        matrix `EV'[1,`j']=`val'
    }
    ereturn matrix sa_bgt=`B'
    ereturn matrix sa_Vgt=`V'
    ereturn matrix sa_Vif=`VI'
    ereturn matrix sa_cells=`C'
    ereturn matrix sa_loading=`L'
    ereturn matrix sa_events=`EV'
    ereturn scalar sa_vce_factor=scalar(`factor')
    ereturn scalar sa_covariance_eigenratio=scalar(`pdratio')
    ereturn scalar sa_N_units=`nu'
    ereturn scalar sa_N_time=`nt'
    ereturn scalar sa_N_groups=`ng'
    ereturn scalar sa_level=`level'
    ereturn local sa_ivar `ivar'
    ereturn local sa_tvar `tvar'
    ereturn local sa_gvar `gvar'
    ereturn local sa_outcome `y'
    ereturn local sa_source_file `"`sourcefile'"'
    ereturn local sa_source_version "0.1"
    ereturn local sa_source_date "24jan2022"
    ereturn local sa_contract "0.6.0-sa1"
    ereturn local sa_weight_scope "fixed"
    ereturn local cmd "jwselect_sa"
    ereturn local version "0.7.0"
    // Empty e(predict)/e(estat_cmd) would invoke Stata's generic defaults on
    // regression columns whose generated variables no longer exist.
    ereturn local predict "jwselect_p"
    ereturn local estat_cmd "jwselect_estat"
    ereturn local marginsnotok "_ALL"
end

program define _jwselect_sa_display
    version 16.1
    di as txt _n "Sun-Abraham cohort-time estimates prepared for jwselect"
    di as txt "Independent panel units: " as res %9.0g e(sa_N_units) ///
        as txt "    Treatment cohorts: " as res %5.0g e(sa_N_groups)
    di as txt "Balanced panel; never-treated controls; full unbinned event indicators; event -1 omitted."
    di as txt "Full unit-clustered regression covariance retained. Aggregation in jwselect uses fixed weights."
    di as txt "Run jwselect, select(wbic), or store these source estimates before comparing specifications."
end
