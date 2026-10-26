*! jwselect 0.5.0: installed-package smoke check, 24oct2026.
*! Run in a FRESH Stata 17+ session from a separate writable folder.
*! Clears data and estimates. Save your work first.
*! Example (replace both paths):
*! cd "C:/work/install_check"
*! do "C:/downloads/jwselect-0.5.0/run_installed_check.do" "C:/downloads/jwselect-0.5.0"
*! The argument is the absolute path of the extracted SOURCE package.
*! This script never changes adopath, installs packages, or reads source fixtures.
version 17.0
args sourcefolder
if `"`sourcefolder'"' == "" {
    display as error "Supply the absolute path of the extracted source package."
    display as error "Run this do-file from a separate writable working directory."
    exit 198
}
capture confirm file `"`sourcefolder'/jwselect.pkg"'
if _rc {
    display as error "The supplied source folder does not contain jwselect.pkg."
    exit 601
}
clear all
set more off
set linesize 120

// Keep all substantive checks inside a captured program so the log closes
// on either success or failure and no caught failure can produce PASS.
program define JW_INSTALL_CHECK
    version 17.0
    syntax , SOURCE(string)
    local source_norm = lower(subinstr(`"`source'"', char(92), "/", .))
    while strlen(`"`source_norm'"')>1 & substr(`"`source_norm'"',-1,1)=="/" {
        local source_norm = substr(`"`source_norm'"',1,strlen(`"`source_norm'"')-1)
    }
    if !(substr(`"`source_norm'"',1,1)=="/" | regexm(`"`source_norm'"',"^[a-z]:/")) {
        display as error "Source path must be absolute."
        exit 198
    }
    if strpos(`"`source_norm'/"',"/../") | strpos(`"`source_norm'/"',"/./") {
        display as error "Use a source path without . or .. directory components."
        exit 198
    }
    capture confirm file "jwselect.ado"
    if !_rc {
        display as error "Current directory contains jwselect.ado; change to a clean working directory."
        exit 459
    }

    display as text "Working directory: `c(pwd)'"
    display as text `"Excluded source folder: `source'"'
    which jwselect, all
    // These are every executable/help file named by jwselect.pkg. Their
    // resolved locations must be outside the supplied source directory.
    foreach required in jwselect.ado _jwselect_jwdid.ado _jwselect_csdid.ado _jwselect_display.ado ///
        _jwselect_mata.do jwselect_estat.ado jwselect_p.ado jwselect.sthlp ///
        _jwselect_load.ado _jwselect_jobguard.ado _jwselect_lockrelease.ado ///
        jwselect_worker.ado jwselect_collect.ado jwselect_worker.sthlp jwselect_collect.sthlp {
        quietly findfile `required'
        local resolved `"`r(fn)'"'
        local resolved_norm = lower(subinstr(`"`resolved'"',char(92),"/",.))
        display as text `"Resolved `required': `resolved'"'
        if !(substr(`"`resolved_norm'"',1,1)=="/" | regexm(`"`resolved_norm'"',"^[a-z]:/")) {
            display as error "Resolved path is not absolute; check adopath and use an installed copy."
            exit 459
        }
        if strpos(`"`resolved_norm'/"',"/../") | strpos(`"`resolved_norm'/"',"/./") {
            display as error "Resolved path contains . or .. components; check adopath."
            exit 459
        }
        if substr(`"`resolved_norm'"',1,strlen(`"`source_norm'"')+1)==`"`source_norm'/"' {
            display as error "A required file resolves to the source folder instead of the installation."
            exit 459
        }
    }

    // Same small known numerical example as the release regression suite.
    // V is a supplied sampling covariance, with structural zero baseline
    // rows. This checks installation and compatible Wald/MD J arithmetic;
    // it does not estimate a covariance matrix from empirical data.
    matrix JIC_B = (.5, 0, .8, .5, 0, 1)
    matrix JIC_C = (3,1 \ 3,2 \ 3,3 \ 4,2 \ 4,3 \ 4,4)
    matrix JIC_V = J(6,6,0)
    matrix JIC_V[1,1] = (.04,0,.0084,.036,0,.0224)
    matrix JIC_V[3,1] = (.0084,0,.05184,.0068,0,.0144)
    matrix JIC_V[4,1] = (.036,0,.0068,.04,0,.0232)
    matrix JIC_V[6,1] = (.0224,0,.0144,.0232,0,.06376)
    matrix JIC_EXPECTED_V = (.05,.01,.03 \ .01,.05,.03 \ .03,.03,.03)
    set obs 1
    generate byte installed_check_marker = 37

    jwselect, bmat(JIC_B) vmat(JIC_V) cells(JIC_C) nunits(100) ///
        select(wbic) wbicscope(joint) candidates(all) maxjoint(4)
    assert "`e(version)'"=="0.5.0"
    assert "`e(build)'"=="0.5.0-adapters1"
    assert "`e(wbic_scope)'"=="joint"
    assert e(G)==2 & e(P)==2 & e(K)==2
    mata: assert(max(abs(st_matrix("e(b)")-(.7,.7,.7)))<1e-10)
    mata: assert(max(abs(st_matrix("e(V)")-st_matrix("JIC_EXPECTED_V")))<1e-10)
    mata: assert(all(st_matrix("e(selection)"):==1))
    estat pretrend
    assert r(df)==2
    assert abs(r(chi2)-125/19)<1e-10
    assert abs(r(J)-125/19)<1e-10
    assert abs(r(p)-chi2tail(2,125/19))<1e-10
    estat pretrend, selected
    assert r(df)==2 & r(postselection)==1
    assert abs(r(chi2)-125/19)<1e-10

    jwselect, bmat(JIC_B) vmat(JIC_V) cells(JIC_C) nunits(100) select(manual)
    mata: assert(max(abs(st_matrix("e(b)")-(.8,1,.9)))<1e-10)
    mata: assert(all(st_matrix("e(selection)"):==0))
    estat pretrend, selected
    assert r(df)==0 & r(chi2)==0 & r(J)==0 & missing(r(p))
    assert _N==1 & installed_check_marker==37
end

local stamp = strtrim(string(clock("`c(current_date)' `c(current_time)'", "DMYhms"), "%18.0f"))
local logfile "jwselect_installed_check_`stamp'.log"
capture log close jwinstall
log using `"`logfile'"', text name(jwinstall)
display as text "Installed-package check: jwselect 0.5.0-adapters1"
display as text "Started: `c(current_date)' `c(current_time)'"
capture noisily JW_INSTALL_CHECK, source(`"`sourcefolder'"')
local check_rc = _rc
if `check_rc' {
    display as error "jwselect 0.5.0 installed-package check: FAIL (return code `check_rc')"
}
else display as result "jwselect 0.5.0 installed-package check: PASS"
display as text "Log: `c(pwd)'/`logfile'"
log close jwinstall
capture program drop JW_INSTALL_CHECK
if `check_rc' exit `check_rc'
