# sdid
Synthetic Difference-in-Differences: Stata package (Beta version)

This Stata package implements the synthetic difference-in-differences estimation procedure, along with a range of inference procedures, following Arkhangelsky et al., (2021).  Arkhangelsky et al. provide a code implementation in R, with accompanying materials here: [synthdid](https://synth-inference.github.io/synthdid/). 
Here we provide a native Stata implementation, principally written in Mata.  This package is currently under active development.

## Inputs
+ Y: Outcome variable (numeric)
+ S: Unit variable (numeric)
+ T: Time variable (numeric)
+ D: Dummy of treatement (numeric)

## Syntax
```
sdid Y S T D, vce(type) seed() reps()
```
+ vce(): bootstrap and placebo (trial version) standard error. 
+ see(): seed define for pseudo-random numbers.
+ reps(): repetitions for bootstrap and placebo se.

**_NOTE:_**  reps() option are required for bootstrap and placebo standard error.


### References
Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. Synthetic Difference in Differences, American Economic Review, December 2021.
