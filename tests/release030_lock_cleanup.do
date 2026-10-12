*! Real filesystem checks for bounded, nonrecursive owned-lock cleanup.
*! Run from the package root. Does not require jwselect Mata to be loaded.
version 17.0
adopath ++ "."
preserve
local saved_rng `"`c(rngstate)'"'
mata: st_local("loaded_before", strofreal(findexternal("jw_version()")!=NULL))
clear
set obs 3
generate long lock_test_sentinel = 71+_n
char _dta[lock_test_note] "preserve helper caller"
tempfile stem
local root `"`stem'_lock test"'
mkdir `"`root'"'
mkdir `"`root'/empty.lock"'
mkdir `"`root'/sibling.lock"'

_jwselect_lockrelease, path(`"`root'/empty.lock"')
assert r(removed)==1 & r(attempts)==1 & r(last_rc)==0 & r(reappeared)==0
mata: assert(!direxists(st_local("root")+"/empty.lock"))
mata: assert(direxists(st_local("root")+"/sibling.lock"))
assert _N==3 & lock_test_sentinel==71+_n
assert `"`: char _dta[lock_test_note]'"'=="preserve helper caller"
assert `"`c(rngstate)'"'==`"`saved_rng'"'

* An already absent owned path is a successful no-op.
_jwselect_lockrelease, path(`"`root'/empty.lock"')
assert r(removed)==1 & r(attempts)==0 & r(last_rc)==0 & r(reappeared)==0

* A real nonempty directory exercises all five failures without fault hooks.
* Its payload and neighboring directory must survive every attempt.
mkdir `"`root'/nonempty.lock"'
tempname payload
file open `payload' using `"`root'/nonempty.lock/keep.txt"', write text
file write `payload' "This payload must never be deleted by lock cleanup." _n
file close `payload'
quietly checksum `"`root'/nonempty.lock/keep.txt"'
local payload_crc = r(checksum)
local payload_bytes = r(filelen)
_jwselect_lockrelease, path(`"`root'/nonempty.lock"') pause(0)
assert r(removed)==0 & r(attempts)==5 & r(last_rc)!=0 & r(reappeared)==0
mata: assert(direxists(st_local("root")+"/nonempty.lock"))
mata: assert(direxists(st_local("root")+"/sibling.lock"))
quietly checksum `"`root'/nonempty.lock/keep.txt"'
assert r(checksum)==`payload_crc' & r(filelen)==`payload_bytes'
assert _N==3 & lock_test_sentinel==71+_n
assert `"`: char _dta[lock_test_note]'"'=="preserve helper caller"
assert `"`c(rngstate)'"'==`"`saved_rng'"'

* Bound invalid options before touching the supplied directory.
capture noisily _jwselect_lockrelease, path(`"`root'/sibling.lock"') maxattempts(0)
assert _rc==198
mata: assert(direxists(st_local("root")+"/sibling.lock"))

* Only this test removes the known test payload; the helper never does.
erase `"`root'/nonempty.lock/keep.txt"'
_jwselect_lockrelease, path(`"`root'/nonempty.lock"')
assert r(removed)==1 & r(attempts)==1 & r(last_rc)==0
_jwselect_lockrelease, path(`"`root'/sibling.lock"')
assert r(removed)==1 & r(attempts)==1 & r(last_rc)==0
rmdir `"`root'"'
mata: assert((findexternal("jw_version()")!=NULL)==strtoreal(st_local("loaded_before")))
assert `"`c(rngstate)'"'==`"`saved_rng'"'
restore
display as result "Lock cleanup: empty/absent paths, bounded failures, payload/data/RNG preservation: PASS"
