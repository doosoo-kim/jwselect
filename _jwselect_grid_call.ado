*! jwselect 0.8.8 03nov2026
program define _jwselect_grid_call
    version 16.1
    args h vh cross actions masks sigma penalty control result scorecov
    tempfile gridinput gridoutput
    tempname fh
    mata: jw_grid_write_input(st_local("gridinput"),tokens(st_local("h")+" "+st_local("vh")+" "+st_local("cross")+" "+st_local("actions")+" "+st_local("masks")+" "+st_local("sigma")+" "+st_local("penalty")+" "+st_local("control")+" "+st_local("scorecov")))
    local platform linux
    if "`c(os)'"=="Windows" local platform w64
    else if strpos(lower("`c(os)'"),"mac") local platform macos
    capture findfile _jwselect_bb_`platform'.plugin
    if _rc exit 497
    local module `"`r(fn)'"'
    capture program _jwselect_bb_088, plugin using(`"`module'"')
    if _rc!=0 & _rc!=110 exit 497
    plugin call _jwselect_bb_088, `"`gridinput'"' `"`gridoutput'"'
    capture confirm file `"`gridoutput'"'
    if _rc {
        di as err "The native BB solver returned no result."
        exit 459
    }
    file open `fh' using `"`gridoutput'"', read text
    file read `fh' reply
    file close `fh'
    gettoken status reply : reply
    if "`status'"=="UNAVAILABLE" exit 499
    if "`status'"=="UNSUPPORTED" exit 498
    if inlist("`status'","BUDGET","RESOURCE") exit 497
    if "`status'"!="OK" {
        di as err "Adaptive BB worker error`reply'"
        exit 459
    }
    mata: st_matrix(st_local("result"),strtoreal(tokens(st_local("reply"))))
end
