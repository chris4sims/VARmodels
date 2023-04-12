# VARmodels
Estimation, impulse responses with error bands, for reduced form VAR's and for structural VAR's identified through heteroskedasticity.
This is an R source package, so should be installable as a local R package on your computer.  The package is well tested on work flows close to that given in the sequence of examples in inst/examples.R (or, after installation, in the main package directory).  But the package has features meant to make it useful with panel data that are not well tested, so may not work at all.
## Pitfalls
Though most have reasonable defaults, there are a lot of parameters in the program arguments.  Setting them badly can lead to slow convergence and/or nonsense results.  `OwnLagMeans` defaults to values appropriate for persistent time series and needs to be changed when some or all series are not persistent (like asset excess returns).  The `sig` parameter vector needs to match the order of magnitude of forecast error standard deviations; if it is off by a factor of 10 or more, results may be nonsense.  A common symptom of problems with these settings is apparently rapidly explosive impulse responses.
## Dependencies
The `examples.R` file invokes my **optimize** package, though other optimizers could do the job.  It also invokes the `coda` package, which is helpful in dianostics for MCMC results.
## Improvements
Suggestions for bug fixes, efficiency improvements, extensions, better documentation, etc. are welcome.  Send a pull request or open an issue.
