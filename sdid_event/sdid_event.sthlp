{smcl}
{viewerjumpto "Syntax" "yatchew_test##syntax"}{...}
{viewerjumpto "Description" "yatchew_test##description"}{...}
{viewerjumpto "Vignettes" "yatchew_test##vignettes"}{...}
{viewerjumpto "Options" "yatchew_test##options"}{...}
{viewerjumpto "Examples" "yatchew_test##examples"}{...}
{viewerjumpto "Saved results" "yatchew_test##saved_results"}{...}

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
{cmd:disag}
{cmd:vce(}{it:string}{cmd:)}
{cmd:brep(}{it:integer} 50{cmd:)]}
{cmd:method(}{it:string}{cmd:)}
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
{cmd:brep()}: number of bootstrap replications (default = 50).
{p_end}

{p 4 4}
{cmd:method()}: estimation method. Allowed arguments:
{cmd:sdid} (default) for Synthetic DiD, {cmd:did} for
traditional DiD and {cmd:sc} for Synthetic Control.
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
{phang2}{stata sdid_event womparl country year quota, vce(bootstrap)}{p_end}
{phang2}{stata mat res = e(H)}{p_end}
{phang2}{stata svmat res}{p_end}
{phang2}{stata gen id = _n - 1 if !missing(res1)}{p_end}
{phang2}{stata replace res1 = 0 in 1}{p_end}
{phang2}{stata replace res3 = 0 in 1}{p_end}
{phang2}{stata replace res4 = 0 in 1}{p_end}
{phang2}{stata twoway  (rarea res3 res4 id, lc(gs10) fc(gs11)) (scatter res1 id, mc(black)) , legend(off) title(sdid_event) xtitle(Relative time to treatment change) ytitle(Women in Parliament) yline(0, lc(red) lp(-))}{p_end}

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

{marker authors}{...}
{title:Authors}

{cmd:sdid_event}
{p 4 4}
Diego Ciccia, Sciences Po. 
{browse "mailto:diego.ciccia@sciencespo.fr":diego.ciccia@sciencespo.fr}
{p_end}

{cmd:sdid}
{p 4 4}
Damian Clarke, Universidad de Chile.
{browse "mailto:dclarke@fen.uchile.cl":dclarke@fen.uchile.cl}
{p_end}

{p 4 4}
Daniel Pailañir, Universidad de Chile.
{browse "mailto:dpailanir@fen.uchile.cl":dpailanir@fen.uchile.cl}
{p_end}




