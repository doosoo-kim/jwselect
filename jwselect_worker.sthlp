{smcl}
{* *! jwselect 0.4.0 23oct2026}{...}
{title:jwselect_worker -- Run prepared focused-inference search jobs}

{p 4 4 2}
Requires Stata 17. Prepare once using {help jwselect} with
{cmd:select(focused) inference(none) prepare("newdir")}. Each job is one
output, one of four objectives, and one starting point. Eight starts give
32 jobs per active output; these need not run in 32 simultaneous sessions.

{title:Syntax}

{p 8 12 2}
{cmd:jwselect_worker using} {it:"directory"}{cmd:,}
{opt worker(#)} {opt workers(#)}
[{opt outputs(numlist)} {opt objectives(numlist)} {opt starts(numlist)}
{opt jobids(numlist)}
{opt cachemb(#)} {opt progress}]

{title:Options}

{p 4 4 2}
{cmd:worker()} is the worker index, from 1 through {cmd:workers()}.
The latter is the number of workers in the requested partition. Job assignment
is deterministic: job ID minus one modulo the worker count equals worker index
minus one. Run each index in a separate Stata session to execute concurrently.
The command itself does not launch Stata processes or submit scheduler jobs.

{p 4 4 2}
{cmd:outputs()}, {cmd:objectives()}, {cmd:starts()}, and {cmd:jobids()} further
restrict the assigned jobs. {cmd:jobids()} lists original global job IDs in
this preparation. All filters are combined with AND and with the modulo-based
worker partition; filtering never renumbers jobs. Output numbers refer to the prepared output map, not cohort
dates. Objectives are 1 lower quantile, 2 upper quantile, 3 left-tail probability,
and 4 right-tail probability. Starts are numbered in the prepared starting-point
matrix. Filters select existing jobs and do not change search settings.
A filtered worker must still satisfy the worker-partition assignment.

{p 4 4 2}
{cmd:cachemb(64)} is the default search-projection cache allowance in decimal
MB, per active calculation. It is not total process memory. The engine caches
as many initial candidate-projection columns as fit and streams the rest in
original candidate order. There is no all-or-nothing cache threshold.
{cmd:cachemb(0)} forces streaming for every candidate. Use the same cache setting when checking exact serial/divided
agreement. Caching and streaming implement the same formulas, but floating-point
rounding near ties can change a search path; assess differences with tolerances
and numerical diagnostics. Concurrent caches require additional RAM.

{p 4 4 2}
{cmd:progress} reports starts and every 50 objective evaluations. A single
large-family evaluation can still take appreciable time. The frozen candidate families, joint
covariance, common search draws, independent validation draws, starts, and budgets
are retained. The caller's RNG state and posted estimation results are preserved.

{title:Restart, locking, and consistency}

{p 4 4 2}
The directory must be a complete preparation from the same release. Configuration
checks reject incompatible or altered inputs. Each assigned job uses a directory
lock {cmd:job_<id>.lock}. Overlapping assignments are refused; disjoint filtered
workers can coexist. Do not delete a lock while its worker is running. If a worker
crashed, establish that its process has stopped before removing a stale lock.

{p 4 4 2}
Cleanup retries a failed empty-directory removal up to five times, waiting
500 milliseconds between failures, and verifies that the directory is gone.
It never deletes directory contents. If removal reports success but the path
is present again, cleanup stops rather than risk deleting another process's
new lock. Persistent cleanup failures print the exact paths. Successful search
computations remain successful and their saved results are retained, but
collection still refuses remaining job locks. Inspect them only after all
processes using the batch have stopped.

{p 4 4 2}
A saved result has a {cmd:job_<id>.bin.done} completion marker. Completed compatible
jobs are reused on restart. An interrupted job without a completion marker can be
recomputed. Do not edit results or copy result files between job IDs/directories.
Changing worker partitioning between nonoverlapping runs is allowed once the old
workers have stopped. Every required job must eventually finish for collection.

{p 4 4 2}
The evaluation budget is fixed by preparation, per objective across its starts.
For {cmd:optstarts(8) optmaxeval(1600)}, each start has 200 evaluations. Running
more sessions does not increase that budget or certify convergence. More
generally, integer quotient budgets are assigned to every start and the
remainder is allocated one evaluation at a time to the first starts. Finer
start-level calls can repeat cache construction. Batch multiple assigned IDs
in {cmd:jobids()} to reuse the loaded context/cache within a target. The focused
application already groups each session's assigned target IDs, including under
start-level division; it uses fixed round-robin assignments, not a dynamic queue. Budget
exhaustion returns the best searched point; mesh convergence is not a proof
of global extrema.

{title:Examples}

{p 4 4 2}
Open four sessions. Use a different worker index in each:

{phang2}{cmd:jwselect_worker using "C:/work/run01", worker(1) workers(4) cachemb(64) progress}
{phang2}{cmd:jwselect_worker using "C:/work/run01", worker(2) workers(4) cachemb(64) progress}
{phang2}{cmd:jwselect_worker using "C:/work/run01", worker(3) workers(4) cachemb(64) progress}
{phang2}{cmd:jwselect_worker using "C:/work/run01", worker(4) workers(4) cachemb(64) progress}

{p 4 4 2}
Alternatively run one output/objective/start job, without other overlapping workers:

{phang2}{cmd:jwselect_worker using "C:/work/run01", worker(1) workers(1) outputs(1) objectives(2) starts(3) progress}

{p 4 4 2}
Batch example (the IDs refer to the original prepared index):

{phang2}{cmd:jwselect_worker using "C:/work/run01", worker(1) workers(1) jobids(1 2 9 10) cachemb(64) progress}

{p 4 4 2}
After every required job is complete and workers have released their locks:

{phang2}{cmd:jwselect_collect using "C:/work/run01", progress}

{title:Saved results}

{p 4 4 2}
{cmd:r(assigned)} is the count of assigned jobs, including previously completed
jobs that were reusable. {cmd:r(worker)}, {cmd:r(workers)}, and
{cmd:r(directory)} record the invocation. Final effects and diagnostics are
returned by {help jwselect_collect}.

{p 4 4 2}
{cmd:r(lock_failures)} counts assigned job locks that could not be released.
Zero means cleanup observed every owned lock absent; a positive value reports
cleanup failure after successful search computations, not numerical failure.
If a search itself fails, its original return code is preserved through cleanup.

{title:Validation}

{p 4 4 2}
Run {cmd:do run_release_tests.do} from the package root in a fresh Stata 17
session. The test clears data and estimates; save your work first.
Source and independent numerical checks do not substitute for executing
the Stata tests.

{title:Also see}

{p 4 4 2}
{help jwselect}, {help jwselect_collect}, {cmd:docs/jwselect_quick_user_guide.pdf}.
