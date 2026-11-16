*! jwselect 0.9.1 16nov2026
program define _jwselect_load, rclass
    version 16.1
    quietly findfile _jwselect_mata.do
    local loader `"`r(fn)'"'
    capture mata: assert(jw_version()=="0.9.1" & jw_build()=="0.9.1-20261116")
    if _rc {
        capture mata: assert(findexternal("jw_version()")!=NULL)
        if !_rc {
            di as error "A different jwselect build is loaded. Start a fresh Stata session."
            exit 499
        }
        capture noisily do `"`loader'"'
        if _rc exit _rc
        mata: assert(jw_version()=="0.9.1" & jw_build()=="0.9.1-20261116")
    }
    return local loader `"`loader'"'
end
