# sdid_event

A Stata module to implement event study analysis with `sdid`.

## Overview

`sdid_event` computes the event-study version of the Synthetic Difference-in-Differences (SDID) estimators from Arkhangelsky et al. (2021). As the name suggests, this program is an extension of `sdid`. As a result, all the conventions and requirement for the implementation of `sdid` also apply for `sdid_event`. `sdid_event` can also be used in staggered adoption designs, with differential timing of treatment adoption. Also, the command supports both the **bootstrap** and **placebo** vce() options for inference.

The derivation of the estimators computed by sdid_event can be found in the [companion paper](https://arxiv.org/abs/2407.09565).
The user can also request cohort-specific aggregated and event study estimates via the **disag** option.

This package depends on `sdid` and `unique`, which can be both installed from SSC.

## Setup

```stata
net install sdid_event, from("https://raw.githubusercontent.com/DiegoCiccia/sdid/main/sdid_event") replace
```

## Syntax

```stata
sdid_event Y G T D [if] [in] [, effects(integer 0) disag vce(string) brep(integer 50)]
```

where:
+ **Y** is the outcome variable.
+ **G** is the unit/group variable.
+ **T** is the time variable.
+ **D** is the treatment variable.

As in `sdid`, the dataset has to be a balanced panel and **D** has to be a binary and absorbing treatment, meaning that the treated units cannot revert their treatment status.

### Options
+ **effects**: number of event study effects to be reported.
+ **disag**: reports estimates of the cohort-specific aggregated and event study estimators.
+ **vce(** off | bootstrap | placebo **)**: selects method for bootstrap inference. With **off**, the program reports only the point estimates, while **bootstrap** and **placebo** correspond to Algorithms 2 and 4 in Clarke et al. (2023).
+ **brep()**: number of bootstrap replications (default = 50).
+ **method(**sdid | did | sc **)**: selects estimation method as in **sdid**.

## Output

+ **e(H)**: matrix with console output.
+ **e(H_c)**: matrix with **disag** option output.
+ **e(b)** and **e(V)**: conventional point estimate vector and variance matrix to allow for integration with **estout**.

## Example

DGP with time-varying treatment effect:

```stata
clear
local GG = 19
local TT = 20
set seed 0
set obs `=`GG' * `TT''

gen G = mod(_n-1,`GG') + 1
bys G: gen T = _n
gen D = T > mod(G, 4) + 1 & G >= `GG'/4
gen Y = uniform() * (1 + D + 10*D*T)

sdid_event Y G T D
```

Generating a graph
```stata
mat res = e(H)
svmat res
gen id = _n - 1 if !missing(res1)

* Turning the ATT line into the reference period
foreach v of varlist res* {
    replace `v' = 0 in 1
}

twoway (line res1 id, lc(black)) (rcap res3 res4 id, lc(black)) (scatter res1 id, mc(black)) , legend(off) title("sdid_event") xtitle("Relative time to treatment change") 
```

The result should look like this:

![sdid_event](https://github.com/DiegoCiccia/sdid/assets/71022390/08917647-08b6-4a52-a0f5-64170bee45ce)

## References 

Arkhangelsky, D., Athey, S., Hirshberg, D., Imbens, G., Wager, S. (2019) [Synthetic difference in differences](https://www.nber.org/papers/w25532)

Ciccia, D (2024) [A Short Note on Event-Study Synthetic Difference-in-Differences Estimators](https://arxiv.org/abs/2407.09565)

Clarke, D. Pailanir, D. Athey, S., Imbens, G. (2023) [Synthetic difference in differences estimation](https://arxiv.org/abs/2301.11859)
