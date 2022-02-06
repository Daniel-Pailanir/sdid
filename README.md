# sdid
Synthetic Difference in Difference :  Stata package

## Inputs
+ Y: Outcome variable (numeric)
+ S: Unit variable (numeric)
+ T: Time variable (numeric)
+ D: Dummy of treatement (numeric)

## Syntax
sdid Y S T D, [seed() breps()]

+ *seed()* and *breps()* option, are required for bootstrap standard error.

### References
Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. Synthetic Difference in Differences, 2019.
