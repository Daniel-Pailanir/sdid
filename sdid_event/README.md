# sdid_event

A Stata module to implement event study analysis with `sdid`.

# Setup

```stata
net install sdid_event, from("https://raw.githubusercontent.com/DiegoCiccia/sdid/main/sdid_event") replace
```

# Adapting `sdid` to event study analysis

In what follows, I use the notation from Clark et al. (2023) to present the estimation procedure for event-study Synthetic Difference-in-Differences (SDID) estimators.
In a setting with $G$ groups observed over $T$ periods, $N_{tr} < G$ groups receive treatment $D$ starting from period $1 < a \leq T$, henceforth called *cohort* or adoption period.
Units are indexed by $i$, such that the first $N_{co} = N - N_{tr}$ units are the controls.
The treatment $D$ is binary, i.e. $D \in \lbrace 0,1\rbrace$, and it affects some outcome of interest $Y$.
The values of $a$ are collected in $A$, i.e. the adoption date vector.
For the sake of generality, we assume that $|A| > 1$, meaning that groups start receiving the treatment at different periods.
The case with no differential timing can be simply retrieved by considering one cohort at a time.
The number of periods from the the onset of the treatment to end of the panel is denoted as $T^{a}_{post}$ and it is cohort-specific.

The cohort-specific SDID estimator from Arkhangelsky et al. (2019) can be rearranged as follows:
$$
\tau^{sdid}_a = \frac{1}{T^a_{post}} \sum_{t = a}^T \left( \frac{1}{N_{tr}} \sum_{i = N_{co} + 1}^N Y_{i,t} - \sum_{i = 1}^{N_{co}} \omega_i Y_{i,t}\right) -  \sum_{t = 1}^{a-1} \left( \frac{1}{N_{tr}} \sum_{i = N_{co} + 1}^N \lambda_t Y_{i,t} - \sum_{i = 1}^{N_{co}}\omega_i \lambda_t  Y_{i,t}\right)
$$
where $\lambda_t$ and $\omega_i$ are the optimal weights chosen to best approximate the pre-treatment outcome evolution of treated and (synthetic) control units.
$\tau^{sdid}_a$ compares the average outcome difference of treated in the $a$ cohort and controls before and after the onset of the treatment.
In doing so, $\tau^{sdid}_a$ encompasses all the post-treatment periods.
As a result, it is possible to estimate the treatment effect $\ell$ periods after the adoption of the treatment, with $\ell \in \lbrace 1,..., T^a_{post} \rbrace$, via a simple disaggregation of $\tau^{sdid}_a$ into the following event-study estimators:

$$
\tau^{sdid}_{a, \ell} = \frac{1}{N_{tr}} \sum_{i = N_{co} + 1}^N Y_{i,a-1+\ell} - \sum_{i = 1}^{N_{co}} \omega_i Y_{i,a-1+\ell} -  \sum_{t = 1}^{a-1} \left( \frac{1}{N_{tr}} \sum_{i = N_{co} + 1}^N \lambda_t Y_{i,t} - \sum_{i = 1}^{N_{co}}  \omega_i \lambda_t Y_{i,t}\right)
$$

This estimator is very similar to those suggested by Borusyak et al. (2021), Liu et al. (2021) and Gardner (2021), when the design is a canonical DiD (de Chaisemartin and D'Haultfoeuille, 2024). The only difference lies in the presence of unit-time specific weights. Notice that by construction:

$$
\tau^{sdid}_a = \frac{1}{T^a_{post}} \sum_{\ell = 1}^{T^a_{post}} \tau^{sdid}_{a, \ell}
$$











# Example

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



