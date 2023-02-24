{smcl}
{* *! version 2.0.0 February 23, 2024}
{title:Title}

{p 4 4 2}
{cmdab:sdid} {hline 2} Synthetic difference-in-differences estimation, inference, and visualization

{marker syntax}{...}
{title:Syntax}

{p 4 4 2}
{opt sdid} {opt depvar} {opt groupvar} {opt timevar} {opt treatment} {ifin}{cmd:,} {it:vce(vcetype)} [{it:options}]

{synoptset 29 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opth vce(vcetype)}}{it: vcetype} may be {opt bootstrap}, {opt jackknife}, or {opt placebo}.{p_end}
{synopt :{opt covariates}({it:{help varlist:varlist}}, [{it:type}])} Allows for the inclusion of covariates in the calculation of the synthetic counterfactual.
Optional {it:type} can be specified, as either "optimized" (the default) or "projected", which is preferable in certain circumstances. {p_end}
{synopt :{opt seed}({it:#})} set random-number seed to #.{p_end}
{synopt :{opt reps}({it:#})} repetitions for bootstrap and placebo inference.{p_end}
{synopt :{opt method(type)}} Allows for an estimation method to be requested.  {it: type} can be "sdid" (which is estimated by default), "did" (for standard difference-in-differencs) or "sc" (for standard synthetic control).{p_end}
{synopt :{opt graph}} if this option is specified, graphs will be displayed in the style of figure 1 from {help sdid##SDID2021:Arkhangelsky et al. (2021)}{p_end}
{synopt :{opt g1on}} If graphing is requested, this option activates the unit-specific weight graph.{p_end}
{synopt :{opt g1_opt}({it:{help twoway_options:graph options}})} option to modify the appearance of the unit-specific weight graph.{p_end}
{synopt :{opt g2_opt}({it:{help twoway_options:graph options}})} option to modify the appearance of the outcome trend graphs.{p_end}
{synopt :{opt graph_export}({it:string}, {it:{help graph export:type}})} option allowing for generated graphs to be saved to the disk.{p_end}
{synopt :{opt msize(markersizestyle)}}{it: markersizestyle} allows you to modify the size of the marker for graph 1.{p_end}
{synopt :{opt unstandardized}} In the case of "optimized" covariates, by default covariates will be standardized as z-scores,
unless the unstandardized option is specified.{p_end}
{synopt :{opt mattitles}} Requests that weights returned in matrices are accompanied by the name
of the groupvar corresponding to each weight.{p_end}
{pstd}
{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}
{pstd}

{pstd}
 {cmd:sdid} implements the synthetic difference-in-differences estimation procedure, along with a range of inference and graphing procedures as described in {help sdid##SDID2021:Arkhangelsky et al. (2021)}.
 Synthetic difference-in-differences is based on a panel (group by time) set-up, in which certain units are treated and
 remaining units are untreated.  The {cmd:sdid} procedure calculates a treatment effect as the pre- versus post-
 difference-in-difference
 between treated units and synthetic control units, where synthetic control units are chosen as an optimally weighted function
 of untreated units (unit-specific weights) and pre-treatment times (time-specific weights).  The {cmd:sdid} command exactly implements
 the procedures described in Arkhangelsky et al. (2021).  The exact estimation procedure implemented by {cmd:sdid} is described in their Algorithm 1.
{p_end}

{pstd}
Much of Arkhangelsky et al. (2021) focuses on cases with a single time period of adoption, however their Appendix A lays out the
estimation procedure in cases of staggered-adoption designs, where treated units can adopt treatment at different moments of
time, while control units never adopt.  {cmd:sdid} seamlessly estimates treatment effects in cases with both single-time periods of treatment and
multiple-time periods of treatment.  In the latter case, rather than calculating a single unit and time-specific weight vector,
an optimal unit and time-specific weight vector is calculated for each adoption period.  The reported average treatment effect
on the treated (ATT) in the staggered adoption design is the weighted estimand described in Arkhangelsky et al. (2021), Appendix A.
Additionally, adoption-period specific estimates and their standard errors are returned after estimation.
{p_end} 

{pstd}
Inference in {cmd:sdid} is based on bootstrap, jackknife, or placebo procedures.  Each procedure is clustered by {opt groupvar},
and follows the precise algorithms laid out in Arkhangelsky et al. (2021).  Specifically, bootstrap inference follows Algorithm 2,
jackknife inference follows Algorithm 3, and placebo inference follows Algorithm 4.  The suitability of each inference procedure
depends on the precise data structure.  For example, bootstrap and jackknife procedures are not appropriate with single treated
units, while placebo inference requires at 1 more control than treated unit.  Inference procedures are provided as standard for
both single-treatment and staggered adoption designs.  In the case of staggered-adoption designs, resample inference is conducted
over the entire ATT, so provides a valid standard error for the headline treatment effect (under the large sample conditions
laid out in Arkhangelsky et al. (2021).)  In very large databases, bootstrap procedures may be computationally expensive, in
which case jackknife will provide a more (computationally) feasible inference procedure.
{p_end}

{pstd}
{cmd:sdid} additionally allows for the inclusion of covariates in a number of ways, and, if requested, provides graphical
output documenting optimal weights, as well as matched treatment and synthetic control trends underlying the synthetic
difference-in-differences framework.  Details related to covariates and graphical options are described at more length below.
{p_end}

{pstd}
Further discussion of this procedure and its implementation in Stata is available in "Synthetic Difference-in-Differences Estimation" (Clarke et al., 2023).
{p_end}

{marker options}{...}
{title:Options}
{dlgtab:Main}
{phang}
{opt vce(vcetype)} is a required option. This may be either bootstrap, jackknife, or placebo, where in each case inference
proceeds following the specified method.  In the case of bootstrap, this is only permitted if greater than one unit is treated.
In the case of jackknife, this is only permitted if greater than one unit is treated in each treatment period (if multiple
treatment periods are considered).  In the case of placebo, this requires at least one more control than treated unit to allow
for permutations to be constructed.  In each case, inference follows the specific algorithm laid out in Arkhangelsky et al. (2021).

{pstd}
{p_end}
{phang}
{opt covariates}({it:{help varlist:varlist}}, [type]) Covariates should be included as a {help varlist:varlist}, and if specified,
treatment and control units will be adjusted based on covariates in the synthetic difference-in-differences procedure.  Optionally,
type may be specified, which indicates how covariate adjustment will occur.  If the type is indicated as "optimized" (the default)
this will follow the method described in Arkhangelsky et al. (2021), footnote 4, where SDID is applied to the residuals of all units
after regression adjustment.  However, this has been observed to be problematic at times (refer to Kranz, 2021), and is also
sensitive to optimization if covariates have high dispersion.  Thus, an alternative type is implmented ("projected"), which consists
of conducting regression adjustment based on parameters estimated only in untreated groups.
This type follows the procedure proposed by Kranz, 2021 (xsynth in R), and is observed to be more stable in some implementations (and at times, considerably faster).
{cmd:sdid} will run simple checks on the covariates indicated and return an error if covariates are constant, to avoid multicolineality.
However, prior to running {cmd:sdid}, you are encouraged to ensure that covariates are not perfectly
multicolinear with other covariates and state and year fixed effects, in a simple two-way fixed
effect regression.  If perfectly multi-colinear covariates are included {cmd:sdid} will execute
without errors, however where type is "optimized", the procedure may be sensitive to the
inclusion of redundant covariates.

{pstd}
{p_end}
{phang}
{opt seed}({it:#}) seed define for pseudo-random numbers.

{pstd}
{p_end}
{phang}
{opt reps}({it:#}) repetitions for bootstrap and placebo se. Default is 50 repetitions.  Larger values should be preferred where possible.

{pstd}
{p_end}
{phang}
{opt method}({it:type}) this option allows for alternative estimation methods to be performed.  Allowed {it:type}s are "sdid" (synthetic difference-in-differences)
"did (standard difference-in-differences) or "sc" (standard synthetic control).  
If this option is not included, sdid is assumed by default.

{pstd}
{p_end}
{phang}
{opt graph} if this option is specified, graphs will be displayed showing unit and time weights as well as outcome trends as per figure 1 from {help sdid##SDID2021:Arkhangelsky et al. (2021)}.

{pstd}
{p_end}
{phang}
{opt g1on} this option activates the unit-specific weight graph. By default g1 is off, as this graph can take considerable time to generate where a large number of control units are present.

{pstd}
{p_end}
{phang}
{opt g1_opt}({it:{help twoway_options:graph options}}) option to modify the appearance of the unit-specific weight graph.
These options adjust the underlying scatter plot, so should be consistent with twoway scatter plots.

{pstd}
{p_end}
{phang}
{opt g2_opt}({it:{help twoway_options:graph options}}) option to modify the appearance of the outcome trend graphs.
These options adjust the underlying line plot, so should be consistent with twoway line plots.

{pstd}
{p_end}
{phang}
{opt graph_export}({it:string}, {it:{help graph export:type}}) Graphs will be saved as weightsYYYY and trendsYYYY for each of the unit-specific weights and outcome trends respectively, where YYYY refers to each treatment adoption period.
Two graphs will be generated for each treatment adoption period. If this option is specified, type must be specified, which refers to a valid Stata graph {it:{help graph export:type}} (eg ".eps", ".pdf", and so forth).
Optionally, a stub can be specified, in which case this will be prepended to exported graph names.

{pstd}
{p_end}
{phang}
{opt msize}({it:{help markersizestyle:size of marker}}) allows you to modify the size of the marker for graph 1.

{pstd}
{p_end}
{phang}
{opt unstandardized} if controls are included and the "optimized" method is specified, controls will be standardized as Z-scores prior to finding optimal weights. This avoids problems with optimization when control variables have very
high dispersion. If unstandardized is specified, controls will simply be entered in their original units.
This option should be used with care.

{pstd}
{p_end}
{phang}
{opt mattitles} Requests labels to be added to the returned {cmd:e(omega)} weight matrix providing names (in string) for the unit variables which generate the synthetic control group in each case.
If mattitles is not indicated, the returned weight matrix ({cmd:e(omega)}) will store these
weights with a final column providing the numerical ID of units, where this numerical ID is either
taken from the unit variable (if this variable is a numerical format), or arranged in alphabetical
order based on the unit variable, if this variable is in string format.

{pstd}
{p_end}


{title:Stored results}

{synoptset 15 tabbed}{...}

{cmd:sdid} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(ATT)}}ATT {p_end}
{synopt:{cmd:e(se)}}Standard error {p_end}
{synopt:{cmd:e(reps)}}Number of bootstrap/placebo replications {p_end}
{synopt:{cmd:e(N_clust)}}Number of clusters {p_end}


{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}sdid{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(vce)}}vcetype specified in vce(){p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}


{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(tau)}}tau estimator for each adoption time-period, along with its standard error{p_end}
{synopt:{cmd:e(lambda)}}lambda weights (time-specific weights){p_end}
{synopt:{cmd:e(omega)}}omega weights (unit-specific weights){p_end}
{synopt:{cmd:e(adoption)}}adoption times{p_end}
{synopt:{cmd:e(beta)}}beta vector of corresponding to covariates (only returned when the covariates option is indicated){p_end}
{synopt:{cmd:e(series)}}control and treatment series of the graphs (only returned when the graph option is indicated) {p_end}
{synopt:{cmd:e(difference)}}difference between treatment and control series (only returned when the graph option is indicated){p_end}

{pstd}
{p_end}

{marker examples}{...}
{title:Examples}

{pstd}
Load data on quotas, women in parliament and maternal mortality (a balanced panel version) from Bhalotra et al., (2020).

{pstd}
 . {stata webuse set www.damianclarke.net/stata/}
 
{pstd}
 . {stata webuse quota_example.dta, clear}

{pstd}
Run SDID estimator without covariates and bootstrap standard error.

{pstd}
 . {stata sdid womparl country year quota, vce(bootstrap) seed(1213)}
 
{pstd}
Run SDID estimator using covariates in projected way.

{pstd}
 . {stata drop if lngdp==.}

{pstd}
 . {stata sdid womparl country year quota, vce(bootstrap) seed(1213) covariates(lngdp, projected)}

{pstd}
Example with one time adoption and some graphics options. Load data from Abadie et al., (2010).

{pstd}
 . {stata webuse prop99_example.dta, clear}

{pstd}
 . {stata sdid packspercapita state year treated, vce(placebo) seed(1213) graph g1_opt(xtitle("")) g2_opt(ylabel(0(50)150))}
 
{marker references}{...}
{title:References}

{marker SDID2021}{...}
{phang} A. Abadie, A. Diamond and J. Hainmueller. 2010. {browse "https://economics.mit.edu/files/11859":{it:Synthetic Control Methods for Comparative Case Studies: Estimating the Effect of California’s Tobacco Control Program}.} Journal of the American Statistical.
{p_end}

{phang}
D. Arkhangelsky, S. Athey, D. Hirshberg, G. Imbens and S. Wager. 2021. {browse "https://www.aeaweb.org/articles?id=10.1257/aer.20190159":{it:Synthetic Difference in Differences}.} 
American Economic Review.
{p_end}

{phang}
S. Bhalotra, D. Clarke, J. Gomes, and A. Venkataramani. 2020. {browse "https://conference.iza.org/conference_files/Gender_2021/clarke_d24360.pdf": {it:Maternal Mortality and Women's Political Participation}.}
Centre for Economic Policy Research Discussion Paper. 
{p_end}

{phang}
D. Clarke, D Pailañir, S. Athey and G. Imbens. 2023. {browse "https://docs.iza.org/dp15907.pdf": {it:Synthetic Difference-in-Differences Estimation}.}
IZA Discussion Paper. 
{p_end}

{phang}
S. Kranz. 2022. {browse "https://github.com/skranz/xsynthdid/blob/main/paper/synthdid_with_covariates.pdf": {it:Synthetic Difference-in-Differences with Time-Varying Covariates}.}
Working paper. 
{p_end}


{title:Author}
Damian Clarke, Universidad de Chile.
Email {browse "mailto:dclarke@fen.uchile.cl":dclarke@fen.uchile.cl}

Daniel Pailañir, Universidad de Chile.
Email {browse "mailto:dpailanir@fen.uchile.cl":dpailanir@fen.uchile.cl}
