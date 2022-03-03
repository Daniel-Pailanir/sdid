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
sdid Y S T D, adoption(type) vce(method) seed() reps()
```
+ adoption(): _**normal**_ and _**staggered**_ adoption; _**normal**_  refers when all units receive the treatment the same period and _**staggered**_ when the treatment is received in different periods of time. It is recommended and can be very helpful to read the appendix of the paper and the repo staggered treatment adoption in R [staggered_adoption_synthdid](https://github.com/zachporreca/staggered_adoption_synthdid).
+ vce(): bootstrap, jackknife and placebo standard error (trial version of bootstrap is avaible now for _**staggered**_ version). 
+ seed(): seed define for pseudo-random numbers.
+ reps(): repetitions for bootstrap and placebo se.

**_NOTE:_**  reps() option are required for bootstrap and placebo standard error.


### References
Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. Synthetic Difference in Differences, American Economic Review, December 2021.
