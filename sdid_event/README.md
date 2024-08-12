# sdid_event

A Stata module to implement event study analysis with `sdid`.

## Overview

`sdid_event` computes the event-study version of the Synthetic Difference-in-Differences (SDID) estimators from Arkhangelsky et al. (2021). As the name suggests, this program is an extension of `sdid`. As a result, all the conventions and requirement for the implementation of `sdid` also apply for `sdid_event`. `sdid_event` can also be used in staggered adoption designs, with differential timing of treatment adoption. Also, the command supports both the **bootstrap** and **placebo** vce() options for inference.

The derivation of the estimators computed by sdid_event can be found in the [companion paper](https://arxiv.org/abs/2407.09565).
The user can also request cohort-specific aggregated and event study estimates via the **disag** option.

This package depends on `sdid` and `unique`, which can be both installed from SSC.

## Setup

### SSC

```stata
ssc install sdid_event, replace
```

### Github

```stata
net install sdid_event, from("https://raw.githubusercontent.com/DiegoCiccia/sdid/main/sdid_event") replace
```

## Syntax

```stata
sdid_event Y G T D [if] [in] [, effects(integer 0) placebo(integer 0) covariates(varlist) disag vce(string) brep(integer 50)]
```

where:
+ **Y** is the outcome variable.
+ **G** is the unit/group variable.
+ **T** is the time variable.
+ **D** is the treatment variable.

As in `sdid`, the dataset has to be a balanced panel and **D** has to be a binary and absorbing treatment, meaning that the treated units cannot revert their treatment status.

### Options
+ **effects**: number of event study effects to be reported.
+ **placebo**: number of placebo estimates to be computed.
+ **covariates**: adds covariates to the estimation routine. To this end, **sdid_event** implements the *projected* method from **sdid**, whereas the outcome is replaced by the residuals of the outcome variable from a TWFE regression on covariates, in the 
sample of untreated and not-yet-treated units.
+ **disag**: reports estimates of the cohort-specific aggregated and event study estimators.
+ **vce(** off | bootstrap | placebo **)**: selects method for bootstrap inference. With **off**, the program reports only the point estimates, while **bootstrap** and **placebo** correspond to Algorithms 2 and 4 in Clarke et al. (2023).
+ **brep()**: number of bootstrap replications (default = 50).
+ **method(** sdid | did | sc **)**: selects estimation method as in **sdid**.

## Output

+ **e(H)**: matrix with console output.
+ **e(H_c)**: matrix with **disag** option output.
+ **e(b)** and **e(V)**: conventional point estimate vector and variance matrix to allow for integration with **estout**.

## Example

```stata
webuse set www.damianclarke.net/stata/
webuse quota_example.dta, clear
keep if quotaYear==2002 | quotaYear==.
sdid_event womparl country year quota, vce(placebo) brep(50) placebo(all)
mat res = e(H)[2..27,1..5]
svmat res
gen id = _n - 1 if !missing(res1)
replace id = 14 - _n if _n > 14 & !missing(res1)
sort id
twoway (rarea res3 res4 id, lc(gs10) fc(gs11%50)) (scatter res1 id, mc(blue) ms(d)), legend(off) title(sdid_event) xtitle(Relative time to treatment change) ytitle(Women in Parliament) yline(0, lc(red) lp(-)) xline(0, lc(black) lp(solid))
```

The result should look like this:

![sdid_event](https://github.com/DiegoCiccia/sdid/assets/71022390/08917647-08b6-4a52-a0f5-64170bee45ce)

## References 

Arkhangelsky, D., Athey, S., Hirshberg, D., Imbens, G., Wager, S. (2019) [Synthetic difference in differences](https://www.nber.org/papers/w25532)

Ciccia, D (2024) [A Short Note on Event-Study Synthetic Difference-in-Differences Estimators](https://arxiv.org/abs/2407.09565)

Clarke, D. Pailanir, D. Athey, S., Imbens, G. (2023) [Synthetic difference in differences estimation](https://arxiv.org/abs/2301.11859)
