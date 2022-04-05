//Specify URL
webuse set www.damianclarke.net/stata/

*------------------------------------------------------------------------------*
*ONE TIME ADOPTION
*------------------------------------------------------------------------------*
webuse quota_example.dta, clear

*-----------*
*Bootstrap SE
*-----------*
//Keep only one time adoption units (2002)
#delimit ;
drop if country=="Algeria" | country=="Jordan" | country=="Kenya"     |
        country=="Rwanda"  | country=="Samoa"  | country=="Swaziland" | 
        country=="Tanzania";
#delimit cr

sdid womparl country year quota, vce(bootstrap) seed(1234) graph

*with control
preserve
drop if lngdp==.
sdid womparl country year quota, vce(bootstrap) seed(1234) covariates(lngdp)
sdid womparl country year quota, vce(bootstrap) seed(1234) covariates(lngdp, projected)
restore

*----------*
*Placebo SE
*----------*
sdid womparl country year quota, vce(placebo) seed(1234) 

*with control: it is necessary drop the missing observations in control variable
preserve
drop if lngdp==.
sdid womparl country year quota, vce(placebo) seed(1234) covariates(lngdp)
sdid womparl country year quota, vce(placebo) seed(1234) covariates(lngdp, projected)
restore

*-----------*
*Jackknife SE
*-----------*
sdid womparl country year quota, vce(jackknife)

*with control: it is necessary drop the missing observations in control variable
preserve
drop if lngdp==.
sdid womparl country year quota, vce(jackknife) covariates(lngdp)
sdid womparl country year quota, vce(jackknife) covariates(lngdp, projected)
restore

*------------------------------------------------------------------------------*
*STAGGERED ADOPTION
*------------------------------------------------------------------------------*
webuse quota_example.dta, clear

*-----------*
*Bootstrap SE
*-----------*
sdid womparl country year quota, vce(bootstrap) seed(1234) graph

*with control: it is necessary drop the missing observations in control variable
preserve
drop if lngdp==.
sdid womparl country year quota, vce(bootstrap) seed(1234) covariates(lngdp)
sdid womparl country year quota, vce(bootstrap) seed(1234) covariates(lngdp, projected)
restore

*----------*
*Placebo SE
*----------*
sdid womparl country year quota, vce(placebo) seed(1234)

*with control: it is necessary drop the missing observations in control variable
preserve
drop if lngdp==.
sdid womparl country year quota, vce(placebo) seed(1234) covariates(lngdp)
sdid womparl country year quota, vce(placebo) seed(1234) covariates(lngdp, projected)
restore

*-----------*
*Jackknife SE
*-----------*
//Keep only treatment periods with more than 2 treated units (2002-2003)
//This is required for jackknife standard error
#delimit ;
drop if country=="Algeria"   | country=="Kenya"    | country=="Samoa"  | 
        country=="Swaziland" | country=="Tanzania";
#delimit cr

sdid womparl country year quota, vce(jackknife)

*with control: it is necessary drop the missing observations in control variable
preserve
drop if lngdp==.
sdid womparl country year quota, vce(jackknife) covariates(lngdp)
sdid womparl country year quota, vce(jackknife) covariates(lngdp, projected)
restore


