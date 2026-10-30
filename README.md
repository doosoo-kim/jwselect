# jwselect 0.7.0

Preperiod selection and minimum-distance estimation in Stata.  
Doosoo Kim · 30 October 2026 · Build `0.7.0-index1`

## Install

Extract the ZIP, then use the folder containing `jwselect.pkg`:

```stata
net install jwselect, from("C:/path/to/jwselect-0.7.0") replace
which jwselect
help jwselect
```

After updating, save your work and restart Stata. Main estimation requires
Stata 16.1 or later; prepared worker/collector commands require Stata 17.
No Python or LaTeX installation is needed.

## Nonlinear conditional-index blocks

After a supported logit or Poisson `jwdid` fit:

```stata
jwdid y x1 x2, ivar(id) tvar(year) gvar(first_treat) ///
    never method(logit) cluster(id)
estimates store NL_SOURCE
jwselect, restrictions(index) select(wbic)
estat design
estat pretrend
```

Each selected date imposes its complete intercept-and-covariate lead block.
Wald–BIC counts the full restriction rank. The joint ATT–index covariance
supports manual, WBIC and focused selection, including reselection inference.
Restore the stored source before changing the selection specification.

Omitting `restrictions()` preserves average-placebo selection. The index
route uses outcome-scale target-level MD; it does not impose natural ATT
bounds. See the guides for identification, covariance and inference conditions.

## Optional upstream packages

The matrix interface needs no upstream estimator. The other routes need only
their relevant optional dependencies; `jwselect` installs none automatically.
Installation instructions and supported specifications are in the guides.

- Linear and nonlinear DiD: `jwdid` with its dependencies.
- CS-DR: `csdid` 2.x, version 2.0.1 or compatible later 2.x. Older SSC 1.x
  interfaces are unsupported.
- Sun–Abraham: `eventstudyinteract`, `reghdfe`, `ftools`, and `avar`; use the
  included `jwselect_sa` wrapper.

## Guides and examples

- [Quick user guide](docs/jwselect_quick_user_guide.pdf)
- [Full user guide](docs/jwselect_user_guide.pdf)
- [Full nonlinear index blocks](examples/nonlinear_index_demo.do)
- [Nonlinear average-placebo selection](examples/nonlinear_demo.do)
- [CS-DR](examples/csdid_demo.do)
- [Sun–Abraham](examples/sa_demo.do)
- [Matrix input](examples/matrix_demo.do)
- [Joint WBIC and pretrend tests](examples/joint_demo.do)
- [Sun–Abraham Monte Carlo do-file](examples/jwselect_sa_monte_carlo.do)
  and [run instructions](examples/jwselect_sa_monte_carlo_README.md)

Run examples from the extracted package folder. Read their opening comments:
some replace the data and estimation results in memory. The Monte Carlo
do-file creates new results only when run.

## Contents and verification

This clean distribution contains the 19 installable Stata files, installation
metadata, two current PDF guides, seven example do-files, Monte Carlo run
instructions, citation metadata, this README, and SHA-256 checksums.
The example datasets are generated or embedded in the do-files.

Historical outputs, saved simulation data, logs, developer tests and build
files, modular source copies, LaTeX sources, and cache files are omitted.
Guide sections describing developer tests refer to the separate full source
archive. The executable Stata files match the validated 0.7.0 release; the
matrix examples are made self-contained for this distribution.

## Citation and contact

Kim, Doosoo. 2026. *From Pretrend Tests to Preperiod Selection in
Difference-in-Differences*. Working paper. Report software version 0.7.0.

Questions and bug reports: [doosoonet@gmail.com](mailto:doosoonet@gmail.com).
Include the Stata version, command, first error, and a reproducible example.
