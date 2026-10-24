ADDITIONAL TESTS FOR JWSELECT 0.4.0
Build: 0.4.0-joint1 | Test revision: 3 | Validation record: 24 October 2026

STATUS

All four additional suites passed across the author's Stata/MP 17 Windows
runs. Supplied logs verify suites 02, 03, and 04. Final PASS for suite 01 was
reported by the author after correcting the test script. The estimation and
postestimation code was unchanged during these additional runs.

RUN THE COMPLETE ADDITIONAL SUITE

Save your work and use a fresh Stata 17+ session. The driver clears data,
matrices, loaded programs, and stored estimates. It uses generated data only.
From the extracted folder containing jwselect.pkg:

   do run_extra_tests.do

The folder is put first on adopath, so this tests the bundled source files.
It does not install or update packages. The separate run_installed_check.do
checks installed files outside the source directory; see README.md.

Results are written to a new tests_extra/results_<timestamp> folder:
   extra_tests.log   Complete plain-text output and version information.
   summary.csv       PASS/FAIL and return code for each of four suites.

The complete-suite success message is:
   jwselect 0.4.0 additional tests: PASS

Some checks deliberately submit invalid inputs and display captured errors.
A suite must reach its final PASS marker to succeed. Any failed suite makes
the overall run FAIL and returns 459. Missing dependencies stop the run with
499 and are never counted as a successful test. Later independent suites
still run after a suite failure, with each result recorded in the summary.

You can choose a NEW output folder whose parent already exists:
   do run_extra_tests.do "tests_extra/my_next_run"

COVERAGE

01_jwdid_joint.do: 17 actual-fit cases
  Generated balanced panel; method(regress) with unconditional/conditional
  covariance, cohort/common baselines, anticipation 0/1, nested/all families;
  plus one default reghdfe fit with conditional covariance. Compares direct
  postestimation with the extracted-matrix route. Independently checks maps,
  joint WBIC, effects/covariance, and all/selected pretrend tests. Also checks
  source data, estimates, sample indicators, and random-number state.

02_design_invariance.do: 35 cases
  Baselines, anticipation 0/1/2, gaps 1/2/0.25, dates after a common baseline,
  equivalent full-set manual restrictions under rebasing, cell permutation,
  raw normalization, outcome rescaling, weights, post windows, aggregate
  covariance, and invalid designs. No claim is made that baseline changes
  must preserve WBIC selection.

03_numeric_oracle.do: 24 designs / 48 fits
  Two/three cohorts with cross-cohort covariance, 960 joint candidate-score
  comparisons, and 96 pretrend tests. Checks masks, effects, loadings, cell
  and reported covariance for joint and cohort WBIC. Empty through fully
  selected sets are included, as are cases where the two selectors differ.
  Reference J independently minimizes the full covariance-weighted MD
  criterion; it is not copied from the Wald statistic. Synthetic covariance
  matrices test the numerical identity, not a covariance estimator's
  robustness. Actual-fit tests separately exercise the robust covariance
  workflow. No Python is needed to execute the Stata tests.

04_saved_results.do: four scenarios
  Joint WBIC, cohort WBIC, manual selection, and focused selection with
  inference(none). Checks stored/disk-loaded results, replay, every documented
  estat subcommand, the prediction guard, and matrix/data/RNG preservation.

REQUIREMENTS AND INDIVIDUAL RUNS

The full driver requires jwdid, jwdid_estat, reghdfe, hdfe, and ftools already
installed. Validated versions are recorded in the supplied logs; no automatic
dependency changes are made. Suites 02-04 use supplied matrices and need no
jwdid. Individual files do not open logs. From the source root, first run:
   adopath ++ "."

Then, for example:
   do tests_extra/01_jwdid_joint.do
   do tests_extra/02_design_invariance.do
   do tests_extra/03_numeric_oracle.do
   do tests_extra/04_saved_results.do

Suites 01, 02, and 04 replace the dataset with synthetic/sentinel data.
Suite 03 leaves the current dataset intact. For a logged rerun of suite 01:
   do run_extra_retry.do

numeric_oracle_fixture.do is input, not a standalone test. Its developer-only
generator, numeric_oracle_reference.py, requires Python, NumPy, and SciPy.
Reference matrices are filled row by row to respect Stata expression limits.

These checks assess implementation behavior on the tested configurations.
They do not establish universal correctness, confidence-interval coverage,
global optimization, or valid chi-squared inference after selection. The
earlier release driver remains the check for focused worker/collector flows.
