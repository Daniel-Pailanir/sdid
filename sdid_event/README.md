# sdid_event

A Stata module to implement event study analysis with `sdid`.

## Setup

```stata
net install sdid_event, from("https://raw.githubusercontent.com/DiegoCiccia/sdid/main/sdid_event") replace
```



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
```

Basic syntax - all the effects are retrieved:

```stata
sdid_event Y G T D
```

**effects()** - limiting the number of reported estimates:

```stata
sdid_event Y G T D, effects(5)
```

**disag** - reporting cohort-specific dynamic treatment effect estimates:

```stata
sdid_event Y G T D, effects(5) disag
```

## References 

Arkhangelsky, D., Athey, S., Hirshberg, D., Imbens, G., Wager, S. (2019) [Synthetic difference in differences](https://www.nber.org/papers/w25532)

Borusyak, K., Jaravel, X. and Spiess, J. (2021) [Revisiting event study designs: Robust and efficient estimation](https://doi.org/10.1093/restud/rdae007)

Clarke, D. Pailanir, D. Athey, S., Imbens, G. (2023) [Synthetic difference in differences estimation](https://arxiv.org/abs/2301.11859)

de Chaisemartin, C., D'Haultfoeuille, X. (2024) [Difference-in-Differences for Simple and Complex Natural Experiments](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4487202)

Liu, L., Wang, Y., Xu, Y. (2024). [A practical guide to counterfactual estimators for causal inference with time‐series cross‐sectional data](https://onlinelibrary.wiley.com/doi/full/10.1111/ajps.12723)

Gardner, J. (2021) [Two-stage differences in differences](https://arxiv.org/abs/2207.05943)



