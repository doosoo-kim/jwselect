*! jwselect 0.7.0 30oct2026: read-only Sun-Abraham adapter
program define _jwselect_sa, rclass
    version 16.1
    syntax [, VCE(string)]
    if "`vce'"=="" local vce supplied
    if lower("`vce'")!="supplied" {
        di as err "The Sun-Abraham adapter uses vce(supplied), the full unit-clustered regression covariance"
        exit 198
    }
    if "`e(cmd)'"!="jwselect_sa" {
        di as err "Run jwselect_sa first; raw eventstudyinteract/reghdfe results lack the validated fitting metadata"
        exit 301
    }
    if "`e(sa_contract)'"!="0.6.0-sa1" | "`e(sa_source_version)'"!="0.1" | ///
        "`e(sa_source_date)'"!="24jan2022" | "`e(sa_weight_scope)'"!="fixed" {
        di as err "Unrecognized Sun-Abraham fitting contract; refit with the current jwselect_sa wrapper"
        exit 498
    }
    if `"`e(datasignaturevars)'"'=="" {
        di as err "The original signed Sun-Abraham estimation sample is unavailable; refit jwselect_sa"
        exit 498
    }
    capture checkestimationsample
    if _rc {
        di as err "The data have changed since jwselect_sa was fitted; restore the estimation data or refit"
        exit 459
    }
    local ivar `e(sa_ivar)'
    local tvar `e(sa_tvar)'
    local gvar `e(sa_gvar)'
    local y `e(sa_outcome)'
    if "`ivar'"=="" | "`tvar'"=="" | "`gvar'"=="" | "`y'"=="" {
        di as err "Sun-Abraham fitting metadata are incomplete"
        exit 498
    }
    confirm numeric variable `ivar' `tvar' `gvar' `y'
    local nu=e(sa_N_units)
    local nt=e(sa_N_time)
    local ng=e(sa_N_groups)
    local nobs=e(N)
    local df=e(df_r)
    local lev=e(sa_level)
    if missing(`nu',`nt',`ng',`nobs',`df',`lev') | ///
        `nu'<=1 | `nt'<2 | `ng'<1 | `nobs'!=`nu'*`nt' | ///
        `nu'!=floor(`nu') | `nt'!=floor(`nt') | `ng'!=floor(`ng') | ///
        e(N_clust)!=`nu' | `df'<=0 | `lev'<=0 | `lev'>=100 | "`e(vce)'"!="cluster" {
        di as err "Invalid Sun-Abraham sample or covariance metadata"
        exit 498
    }
    quietly count if e(sample)
    if r(N)!=`nobs' {
        di as err "The Sun-Abraham estimation sample is missing or has changed"
        exit 498
    }
    tempname B V C L B0 V0 T bad asym eig pd IX VV EE
    foreach mat in b V sa_bgt sa_Vgt sa_cells sa_loading {
        capture confirm matrix e(`mat')
        if _rc {
            di as err "The Sun-Abraham source does not contain e(`mat')"
            exit 498
        }
    }
    matrix `B'=e(sa_bgt)
    matrix `V'=e(sa_Vgt)
    matrix `C'=e(sa_cells)
    matrix `L'=e(sa_loading)
    matrix `B0'=e(b)
    matrix `V0'=e(V)
    local m=`ng'*`nt'
    local k=colsof(`B0')
    if rowsof(`B')!=1 | colsof(`B')!=`m' | ///
        rowsof(`V')!=`m' | colsof(`V')!=`m' | ///
        rowsof(`C')!=`m' | colsof(`C')!=2 | ///
        rowsof(`L')!=`m' | colsof(`L')!=`k' | ///
        rowsof(`B0')!=1 | rowsof(`V0')!=`k' | colsof(`V0')!=`k' {
        di as err "Sun-Abraham coefficient, covariance, and cell-map dimensions do not conform"
        exit 498
    }
    mata: st_numscalar("`bad'", hasmissing(st_matrix("`B'")) | hasmissing(st_matrix("`V'")) | hasmissing(st_matrix("`C'")) | hasmissing(st_matrix("`L'")) | hasmissing(st_matrix("`B0'")) | hasmissing(st_matrix("`V0'")))
    if scalar(`bad') {
        di as err "Sun-Abraham inputs contain missing coefficients or covariance entries"
        exit 498
    }
    mata: st_numscalar("`bad'", max(abs(st_matrix("`B'")-st_matrix("`B0'")*st_matrix("`L'")'))/max((1,max(abs(st_matrix("`B'"))))))
    if scalar(`bad')>1e-10 {
        di as err "Sun-Abraham cell estimates do not match the stored underlying regression"
        exit 498
    }
    mata: st_numscalar("`bad'", max(abs(st_matrix("`V'")-st_matrix("`L'")*st_matrix("`V0'")*st_matrix("`L'")'))/max((1,max(abs(st_matrix("`V'"))))))
    mata: st_numscalar("`asym'", max(abs(st_matrix("`V'")-st_matrix("`V'")'))/max((1,max(abs(st_matrix("`V'"))))))
    mata: st_numscalar("`eig'", min(symeigenvalues(st_matrix("`V'")))/max((1,max(abs(st_matrix("`V'"))))))
    if scalar(`bad')>1e-10 | scalar(`asym')>1e-10 | scalar(`eig')< -1e-10 {
        di as err "Sun-Abraham inputs do not contain the valid, complete regression joint covariance"
        exit 498
    }
    quietly levelsof `gvar' if e(sample) & `gvar'>0, local(groups)
    quietly levelsof `tvar' if e(sample), local(times)
    if `: word count `groups''!=`ng' | `: word count `times''!=`nt' {
        di as err "Sun-Abraham sample support does not match the saved cohort-time grid"
        exit 498
    }
    local row 0
    foreach gg of local groups {
        foreach tt of local times {
            local ++row
            if `C'[`row',1]!=`gg' | `C'[`row',2]!=`tt' {
                di as err "Sun-Abraham cohort-time map does not match its validated fitting design"
                exit 498
            }
            if `tt'==`gg'-1 {
                mata: st_numscalar("`bad'", max(abs(st_matrix("`V'")[`row',.])))
                if `B'[1,`row']!=0 | scalar(`bad')!=0 {
                    di as err "Only the event -1 baseline can be the structural-zero normalization"
                    exit 498
                }
            }
            else if `V'[`row',`row']<=0 {
                di as err "A supported nonbaseline Sun-Abraham cell is unidentified"
                exit 498
            }
        }
    }
    mata: `IX'=selectindex(st_matrix("`C'")[,2]:!=st_matrix("`C'")[,1]:-1)
    mata: `VV'=st_matrix("`V'")
    mata: `EE'=symeigenvalues(`VV'[`IX',`IX'])
    mata: st_numscalar("`pd'",min(`EE')/max(abs(`EE')))
    mata: mata drop `IX' `VV' `EE'
    if missing(scalar(`pd')) | scalar(`pd')<=1e-10 {
        di as err "The identified Sun-Abraham cell covariance is singular or numerically ill-conditioned"
        exit 498
    }
    matrix `T'=J(`m',8,.)
    matrix `T'[1,1]=`C'
    local crit=invttail(`df',(100-`lev')/200)
    forvalues j=1/`m' {
        matrix `T'[`j',3]=`B'[1,`j']
        matrix `T'[`j',4]=sqrt(`V'[`j',`j'])
        if `T'[`j',4]>0 {
            matrix `T'[`j',5]=`T'[`j',3]/`T'[`j',4]
            matrix `T'[`j',6]=2*ttail(`df',abs(`T'[`j',5]))
            matrix `T'[`j',7]=`T'[`j',3]-`crit'*`T'[`j',4]
            matrix `T'[`j',8]=`T'[`j',3]+`crit'*`T'[`j',4]
        }
    }
    matrix colnames `T'=cohort time estimate se statistic pvalue ll ul
    return matrix bgt=`B'
    return matrix Vgt=`V'
    return matrix cells=`C'
    return matrix upstream_table=`T'
    return scalar N_units=`nu'
    return scalar N_obs=`nobs'
    return scalar N_clusters=`nu'
    return scalar N_covariates=0
    return scalar gap=1
    return scalar anticipation=0
    return scalar df_r=`df'
    return local ivar `ivar'
    return local tvar `tvar'
    return local gvar `gvar'
    return local sourcecmd `"`e(cmdline)'"'
    return local sourcevce "cluster"
    return local sourceversion "`e(sa_source_version)'"
    return local sourceestimator "eventstudyinteract"
    return local model "linear"
    return local method "sunabraham"
    return local restriction_scale "outcome_contrast"
    return local upstream_stat "t"
    return local vce "supplied"
end
