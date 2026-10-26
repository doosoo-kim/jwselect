*! jwselect 0.5.0 23oct2026
program define _jwselect_load, rclass
    version 16.1
    quietly findfile _jwselect_mata.do
    local loader `"`r(fn)'"'
    capture mata: assert(jw_version()=="0.5.0-adapters1")
    if _rc {
        capture mata: assert(findexternal("jw_version()")!=NULL)
        if !_rc {
            di as error "A different jwselect build is loaded. Start a fresh Stata session."
            exit 499
        }
        capture noisily do `"`loader'"'
        if _rc exit _rc
        mata: assert(jw_version()=="0.5.0-adapters1")
    }
    return local loader `"`loader'"'
end
