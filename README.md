# sdid -- Synthetic Difference-in-Differences for Stata

This Stata package implements the synthetic difference-in-differences estimation procedure, along with a range of inference procedures, following Arkhangelsky et al., (2021).  Arkhangelsky et al. provide a code implementation in R, with accompanying materials here: [synthdid](https://synth-inference.github.io/synthdid/). 
Here we provide a native Stata implementation, principally written in Mata.  This package is currently under active development.

## Inputs
+ Y: Outcome variable (numeric)
+ S: Unit variable (numeric)
+ T: Time variable (numeric)
+ D: Dummy of treatement, equal to 1 if units are treated (numeric)

## Syntax
```s
sdid Y S T D, vce(method) seed(#) reps(#) controls(varlist) control_type() graph(on)
```
+ vce(): bootstrap, jackknife and placebo standard error.
+ seed(): seed define for pseudo-random numbers.
+ reps(): repetitions for bootstrap and placebo se.
+ controls(  varlist [, method]): controls included to adjust Y.  A varlist of controls should be included, and optionally an option for the method used to adjust.  This can be R in which case it follows the method proposed by Arkhangelsky et al., or xsynth, in which case it follows the procedure proposed by xsynth in R.  Where method is not specified, R is used as default.
+ graph(): "_on_" for display graph in figure 1 from Arkhangelsky et al.

### References
Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. Synthetic Difference in Differences, American Economic Review, December 2021.
