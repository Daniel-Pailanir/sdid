![Stata](https://img.shields.io/badge/stata-2013-green) ![GitHub Starts](https://img.shields.io/github/stars/Daniel-Pailanir/sdid?style=social) ![GitHub forks](https://img.shields.io/github/forks/Daniel-Pailanir/sdid?style=social) ![issues](https://img.shields.io/github/issues/Daniel-Pailanir/sdid) ![issues close](https://img.shields.io/github/issues-closed/Daniel-Pailanir/sdid) ![GitHub license](https://img.shields.io/github/license/Daniel-Pailanir/sdid)

# sdid -- Synthetic Difference-in-Differences for Stata

This Stata package implements the synthetic difference-in-differences estimation procedure, along with a range of inference and graphing procedures, following Arkhangelsky et al., (2021).  Arkhangelsky et al. provide a code implementation in R, with accompanying materials here: [synthdid](https://synth-inference.github.io/synthdid/). 
Here we provide a native Stata implementation, principally written in Mata.  This package extends the funcionality of the original R package, allowing very simply for estimation in contexts with staggered adoption over multiple treatment periods (as well as in a single adoption period as in the original code).  Some further details can be found in the accompanying working paper [here](https://docs.iza.org/dp15907.pdf).

>[!Tip]
>If you wish to implement Event Study analysis with SDiD, please check out [sdid_event](https://github.com/Daniel-Pailanir/sdid/tree/main/sdid_event) and its [technical note](https://arxiv.org/abs/2407.09565).
>

To install directly into Stata:
```s
ssc install sdid, replace
```
or using ```net install``` command:
```s
net install sdid, from("https://raw.githubusercontent.com/daniel-pailanir/sdid/master") replace
```
## Inputs
+ Y: Outcome variable (numeric)
+ S: Unit variable (numeric or string)
+ T: Time variable (numeric)
+ D: Dummy of treatement, equal to 1 if units are treated, and otherwise 0 (numeric)

## Syntax
```s
sdid Y S T D [if] [in], vce(method) seed(#) reps(#) covariates(varlist [, method])
                        zeta_lambda(real) zeta_omega(real) min_dec(real) max_iter(real)
                        method(methodtype) unstandardized graph_export([stub] , type) mattitles
                        graph g1on g1_opt(string) g2_opt(string) msize() 
```
+ vce(): **bootstrap**, **jackknife** and **placebo**. If you want to omit this procedure use **noinference**.
+ method(): **sdid** for Synthetic DiD, **did** for DiD and **sc** for Synthetic Control.
+ seed(): seed define for pseudo-random numbers.
+ reps(): repetitions for bootstrap and placebo se.
+ covariates(  varlist [, method]): covariates included to adjust Y.  A varlist of covariates should be included, and optionally an option for the method used to adjust.  This can be "optimized" in which case it follows the method proposed by Arkhangelsky et al., or "projected", in which case it follows the procedure proposed by Kranz, 2021 (xsynth in R).  Where method is not specified, optimized is used as default.  Kranz has shown that the projected method is preferable in a number of circumstances.  In this implementation, the projected method is often considerably faster.
+ _not_yet: If covariates are included with the "projected" method, allows for projections to be based off all not yet treated units, rather than units which are never treated.
+ zeta_lambda: Define the regularization term $\zeta$. Default is 1e-6.
+ zeta_omega: Define the regularization term $\zeta$. Default is $(N_{tr} T_{post})^{1/4}$ for SDID estimator and 1e-6 for the rest of estimators.
+ min_dec: Define the MinDec term in $(\text{MinDec}\times \sigma_t)^2$. Is a stopping criteria for our weight estimator. Default is 1e-5.
+ max_iter: Defines the maximum number of iterations. Default is 10,000.
+ graph: If this option is specified, graphs will be displayed showing unit and time weights as well as outcome trends as per figure 1 from Arkhangelsky et al.
+ g1on: this option activates the unit-specific weight graph. By default g1 is off.
+ g1_opt() and g2_opt(): option to modify the appearance of the aforementioned graphs.  g1 refers to the unit-specific weight graph, while g2 refers to the outcome trend groaphs. Options requested have to follow the syntax of Stata's twoway_options.
+ unstandardized: If controls are included and the "optimized" method is specified, controls will be standardized as Z-scores prior to finding optimal weights.  This avoids problems with optimization when control variables have very high dispersion.  If unstandardized is specified, controls will simply be entered in their original units.  This option should be used with care.
+ msize(): allows you to modify the size of the marker for graph 1.
+ graph_export( [stub], type): Generated graphs can be saved to the disk using this option.  Graphs will be saved as weightsYYYY and trendsYYYY for each of the unit-specific weights and outcome trends respectively, where YYYY refers to each treatment adoption period.  Two graphs will be generated for each treatment adoption period.  If this option is specified, type must be specified, which refers to a valid Stata graph type (eg ".eps", ".pdf", and so forth).  Optionally, a stub can be specified, in which case this will be prepended to exported graph names. 
+ mattitles: Requests labels to be added to the returned Omega weight matrix providing names (in string) for the unit variables which generate the synthetic control group in each case.  If mattitles is not indicated, the returned weight matrix (e(omega)) will store these weights with a final column providing the numerical ID of units, where this numerical ID is either taken from the unit variable (if this variable is a numerical format), or arranged in alphabetical order based on the unit variable, if this variable is in string format.

## Running an example
```s
webuse set www.damianclarke.net/stata/
webuse prop99_example.dta, clear

#delimit ;
sdid packspercapita state year treated, vce(placebo) reps(100) seed(123) 
     graph g1on g1_opt(xtitle("") ylabel(-35(5)10) scheme(plotplainblind)) 
     g2_opt(ylabel(0(50)150) xlabel(1970(5)2000) ytitle("Packs per capita") 
            xtitle("") text(125 1995 "ATT = -15.604" " SE = (9.338)") scheme(plotplainblind))
    graph_export(sdid_, .png);
#delimit cr
```

The code returns the following results

```
Placebo replications (100). This may take some time.
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5
..................................................     50
..................................................     100


Synthetic Difference-in-Differences Estimator

-----------------------------------------------------------------------------
packsperca~a |     ATT     Std. Err.     t      P>|t|    [95% Conf. Interval]
-------------+---------------------------------------------------------------
   treatment | -15.60383    9.33752    -1.67    0.095   -33.90504     2.69738
-----------------------------------------------------------------------------
95% CIs and p-values are based on Large-Sample approximations.
Refer to Arkhangelsky et al., (2020) for theoretical derivations.
(file sdid_weights1989.eps written in EPS format)
(file sdid_trends1989.eps written in EPS format)
```

<img src="https://github.com/Daniel-Pailanir/sdid/blob/main/graphs/sdid_trends1989.png" width="600" height="400">
<img src="https://github.com/Daniel-Pailanir/sdid/blob/main/graphs/sdid_weights1989.png" width="600" height="400">

To export results, you can use ```eststo``` and ```esttab```:

```s
*create a uniform variable to use as a control
gen r=runiform()

*run sdid
eststo sdid_1: sdid packspercapita state year treated, vce(placebo) seed(2022)
eststo sdid_2: sdid packspercapita state year treated, vce(placebo) seed(2022) covariates(r, projected)

*create a table
esttab sdid_1 sdid_2, starlevel ("*" 0.10 "**" 0.05 "***" 0.01) b(%-9.3f) se(%-9.3f)
```

The code returns the following results

```s
--------------------------------------------
                      (1)             (2)   
             packsperca~a    packsperca~a   
--------------------------------------------
treated           -15.604*        -15.750*  
                  (7.981)         (8.039)   
--------------------------------------------
N                    1209            1209   
--------------------------------------------
Standard errors in parentheses
* p<0.10, ** p<0.05, *** p<0.01
```

### Acknowledgments
We are grateful to Noah Spencer, Jared Greathouse and Asjad Naqvi for very useful feedback related to this code.  We are also very grateful to many other users who have suggested a range of useful extensions and filed bug reports which have made this code better.


### References

Arkhangelsky, D., Athey, S., Hirshberg, D., Imbens, G., Wager, S. (2019) [Synthetic difference in differences](https://www.nber.org/papers/w25532), American Economic Review, December 2021.

Ciccia, D. (2024) [A Short Note on Event-Study Synthetic Difference-in-Differences Estimators](https://arxiv.org/abs/2407.09565)

Clarke, D. Pailanir, D. Athey, S., Imbens, G. (2023) [Synthetic difference in differences estimation](https://arxiv.org/abs/2301.11859), IZA Discussion Paper, January 2023.
