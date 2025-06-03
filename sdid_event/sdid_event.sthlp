{smcl}
{viewerjumpto "Syntax" "sdid_event##syntax"}{...}
{viewerjumpto "Description" "sdid_event##description"}{...}
{viewerjumpto "Options" "sdid_event##options"}{...}
{viewerjumpto "Examples" "sdid_event##examples"}{...}
{viewerjumpto "Saved results" "sdid_event##saved_results"}{...}

{title:Title}

{p 4 4}
{cmd:sdid_event} {hline 2} Synthetic Difference-in-Differences (SDID) event-study estimators.
{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 4}
{cmd:sdid_event Y G T D [if] [in]}
{cmd:[,}
{cmd:effects(}{it:integer} 0{cmd:)}
{cmd:placebo(}{it:integer} 0{cmd:)}
{cmd:covariates(}{it:string}{cmd:)}
{cmd:disag}
{cmd:vce(}{it:string}{cmd:)}
{cmd:brep(}{it:integer} 50{cmd:)}
{cmd:method(}{it:string}{cmd:)}
{cmd:combine(}{it:string}{cmd:)}
{cmd:vcov}
{cmd:sb}
{cmd:boot_ci}
{cmd:_not_yet(}{it:string}{cmd:)}
{cmd:unstandardized]}
{p_end}

{p 4 4}
{cmd:Y} is the outcome variable.
{p_end}

{p 4 4}
{cmd:G} is the unit/group variable.
{p_end}

{p 4 4}
{cmd:T} is the time variable.
{p_end}

{p 4 4}
{cmd:D} is the treatment variable.
{p_end}

{marker description}{...}
{title:Description}

{p 4 4}
{cmd:sdid_event} computes the event-study version of 
the Synthetic Difference-in-Differences (SDID) estimators from 
Arkhangelsky et al. (2021). As the name suggests, this program 
is an extension of {cmd:sdid}. As a result, all the conventions 
and requirement for the implementation of {cmd:sdid} also apply 
for {cmd:sdid_event}. {cmd:sdid_event} can also be used in 
staggered adoption designs, with differential timing of 
treatment adoption. Also, the command supports both the 
{cmd:bootstrap} and {cmd:placebo} vce() options for inference.
{p_end}

{p 4 4}
The derivation of the estimators computed by {cmd:sdid_event} 
can be found in the {browse "https://arxiv.org/abs/2407.09565":companion paper}.
The user can also request cohort-specific aggregated 
and event study estimates via the {cmd:disag} option.
As in {cmd:sdid}, the dataset has to be a balanced panel 
and {cmd:D} has to be a binary and absorbing treatment, 
meaning that the treated units cannot revert their treatment status.
{p_end}

{p 4 4}
This package depends on {cmd:sdid} and {cmd:unique}, 
which can be both installed from SSC.
{p_end}

{marker options}{...}
{title:Options}
{p 4 4}
{cmd:effects()}: number of event study effects to be reported.
By default, all feasible dynamic effects are reported.
{p_end}

{p 4 4}
{cmd:placebo()}: number of placebo estimates to be computed. 
{cmd:placebo(all)} returns all feasible placebo estimates.
{p_end}

{p 4 4}
{cmd:disag}: reports estimates of the cohort-specific aggregated and event study estimators.
{p_end}

{p 4 4}
{cmd:vce()}: selects method for bootstrap inference. 
The allowed arguments are {cmd:off}, {cmd:bootstrap} and {cmd:placebo}.
With {cmd:off}, the program reports only the point estimates,
while {cmd:bootstrap} and {cmd:placebo} correspond to 
Algorithms 2 and 4 in Clarke et al. (2023).
{p_end}

{p 4 4}
{cmd:covariates()}: adds covariates to the estimation routine.
{cmd:sdid_event} now supports all covariate adjustment methods
available in {cmd:sdid}. Covariates should be included as a 
{help varlist:varlist}, and optionally a method may be specified:
{cmd:projected} (default) where the outcome is replaced by the residuals 
of the outcome variable from a TWFE regression on covariates, in the 
sample of untreated and not-yet-treated units, following the procedure 
proposed by Kranz (2022); or {cmd:optimized} which follows the method 
described in Arkhangelsky et al. (2021), footnote 4, where SDID is 
applied to the residuals of all units after regression adjustment.
{p_end}

{p 4 4}
{cmd:_not_yet()}: if covariates are included with the {cmd:projected} method, 
allows for projections to be based off all not-yet-treated units rather 
than only never-treated units. This option is enabled by default for the 
projected method to maintain backward compatibility. Can be explicitly 
disabled with {cmd:_not_yet(off)} to use only never-treated units for projection.
Note: This option always requires parentheses.
{p_end}

{p 4 4}
{cmd:unstandardized}: if covariates are included with the {cmd:optimized} method, 
prevents standardization of covariates as z-scores prior to finding optimal weights. 
By default, covariates are standardized when using the optimized method to avoid 
problems with optimization when control variables have very high dispersion.
{p_end}

{p 4 4}
{cmd:brep()}: number of bootstrap replications (default = 50).
{p_end}

{p 4 4}
{cmd:method()}: estimation method. Allowed arguments:
{cmd:sdid} (default) for Synthetic DiD, {cmd:did} for
traditional DiD and {cmd:sc} for Synthetic Control.
{p_end}

{p 4 4}
{cmd:combine()}: grouping multiple event study 
coefficients under a single estimate. For instance,
with year-group data over 6 years, one could be interested in
comparing differential outcomes in three-year windows 
after the start of the treatment. This can be achieved
by {cmd:combine(1 2 3; 4 5 6)}. This option returns the 
corresponding estimate, plus standard errors, CIs and
number of treated units x post treatment periods in the
requested windows.
{p_end}

{p 4 4}
{cmd:vcov}: returns the variance-covariance matrix
of the requested dynamic effects. If {cmd:placebo()}
is requested, the option also returns the 
variance-covariance matrix of the placebo estimates.
{p_end}

{p 4 4}
{cmd:sb}: returns a matrix with the values of 
the requested estimates across all bootstrap repetitions.
{p_end}

{p 4 4}
{cmd:boot_ci}: by default, 95% CIs are computed 
using a normal approximation for the bootstrap 
distribution of the estimates. With this option on, 
the reported CIs are computed using the empirical 
CDF. 
{p_end}

{marker examples}{...}
{title:Examples}

{p 2 4}
Example 1: Random DGP
{p_end}

{phang2}{stata clear}{p_end}
{phang2}{stata local GG = 19}{p_end}
{phang2}{stata local TT = 20}{p_end}
{phang2}{stata set seed 0}{p_end}
{phang2}{stata set obs `=`GG' * `TT''}{p_end}

{phang2}{stata gen G = mod(_n-1,`GG') + 1}{p_end}
{phang2}{stata gen T = floor((_n-1)/`GG') + 1}{p_end}
{phang2}{stata gen D = T > mod(G, 4) + 1 & G >= `GG'/4}{p_end}
{phang2}{stata gen Y = uniform() * (1 + D + 10*D*T)}{p_end}

{phang2}{stata sdid_event Y G T D}{p_end}

{p 2 4}
Example 2: Bhalotra, Clarke, Gomes 
& Venkataramani (2023) from Section 4.4 
of Clarke et al. (2023)
{p_end}

{phang2}{stata webuse set www.damianclarke.net/stata/}{p_end}
{phang2}{stata webuse quota_example.dta, clear}{p_end}
{phang2}{stata keep if quotaYear==2002 | quotaYear==.}{p_end}
{phang2}{stata sdid_event womparl country year quota, vce(placebo) brep(50) placebo(all)}{p_end}
{phang2}{stata mat res = e(H)[2..27,1..5]}{p_end}
{phang2}{stata svmat res}{p_end}
{phang2}{stata gen id = _n - 1 if !missing(res1)}{p_end}
{phang2}{stata replace id = 14 - _n if _n > 14 & !missing(res1)}{p_end}
{phang2}{stata sort id}{p_end}
{phang2}
{stata twoway (rarea res3 res4 id, lc(gs10) fc(gs11%50)) (scatter res1 id, mc(blue) ms(d)), legend(off) title(sdid_event) xtitle(Relative time to treatment change) ytitle(Women in Parliament) yline(0, lc(red) lp(-)) xline(0, lc(black) lp(solid))}
{p_end}

{p 2 4}
Example 3: Using covariates with projected method (default)
{p_end}

{phang2}{stata webuse set www.damianclarke.net/stata/}{p_end}
{phang2}{stata webuse quota_example.dta, clear}{p_end}
{phang2}{stata keep if quotaYear==2002 | quotaYear==.}{p_end}
{phang2}{stata drop if lngdp==.}{p_end}
{phang2}{stata sdid_event womparl country year quota, covariates(lngdp) vce(bootstrap) brep(50)}{p_end}

{p 2 4}
Example 4: Using covariates with optimized method
{p_end}

{phang2}{stata sdid_event womparl country year quota, covariates(lngdp, optimized) vce(bootstrap) brep(50)}{p_end}

{p 2 4}
Example 5: Using projected method with only never-treated units
{p_end}

{phang2}{stata sdid_event womparl country year quota, covariates(lngdp, projected) _not_yet(off) vce(bootstrap) brep(50)}{p_end}

{marker saved_results}{...}
{title:Saved results}

{p 4 4}
{cmd:e(H)}: matrix with console output.
{p_end}

{p 4 4}
{cmd:e(H_c)}: matrix with {cmd:disag} option output.
{p_end}

{p 4 4}
{cmd:e(b)} and {cmd:e(V)}: conventional point estimate vector and variance matrix to allow for integration with {cmd:estout}.
{p_end}

{marker references}{...}
{title:References}

Arkhangelsky, D., Athey, S., Hirshberg, D., Imbens, G., Wager, S. (2019) {browse "https://www.nber.org/papers/w25532":Synthetic difference in differences}.

Ciccia, D. (2024) {browse "https://arxiv.org/abs/2407.09565":A Short Note on Event-Study Synthetic Difference-in-Differences Estimators}.

Clarke, D. Pailanir, D. Athey, S., Imbens, G. (2023) {browse "https://arxiv.org/abs/2301.11859":Synthetic difference in differences estimation}.

Kranz, S. (2022) {browse "https://github.com/skranz/xsynthdid/blob/main/paper/synthdid_with_covariates.pdf":Synthetic Difference-in-Differences with Time-Varying Covariates}.

{marker authors}{...}
{title:Authors}

{p 4 4}
Diego Ciccia, Northwestern University, Kellogg School of Management. 
{browse "mailto:diego.ciccia@kellogg.northwestern.edu":diego.ciccia@kellogg.northwestern.edu}
{p_end}

{p 4 4}
Damian Clarke, Universidad de Chile.
{browse "mailto:dclarke@fen.uchile.cl":dclarke@fen.uchile.cl}
{p_end}

{p 4 4}
Daniel Pailañir, Universidad de Chile.
{browse "mailto:dpailanir@fen.uchile.cl":dpailanir@fen.uchile.cl}
{p_end}