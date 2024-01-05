*! sdid: Synthetic Difference-in-Differences
*! Version 2.0.0 February 24, 2023
*! Author: Pailañir Daniel, Clarke Damian
*! dpailanir@fen.uchile.cl, dclarke@fen.uchile.cl

/*
Versions
1.0.0 Apr 04, 2022: Original version, staggered adoption [SSC]
1.1.0 May 13, 2022: Correction for 13.0 Mata, bug fix, adding if/in
1.2.0 Jul 12, 2022: Exporting additional details         [SSC]
1.3.0 Jul 13, 2022: Addition of DiD and SC methods
1.4.0 Jan 15, 2023: Standard error for each adoption (Bootstrap only)
2.0.0 Feb 24, 2023: Standard error for each adoption (Bootstrap, placebo, jackknife) [SSC]
*/

program sdid, eclass
version 13.0

#delimit ;
    syntax varlist(min=4 max=4) [if] [in], vce(string) 
    [
    seed(numlist integer >0 max=1)
    reps(integer 50)
    covariates(string asis)
    zeta_lambda(real 1e-6)
    zeta_omega(real 1e-6)
    min_dec(real 1e-5)
    max_iter(real 1e4)
    level(integer 95)
    graph
    g1on
    g1_opt(string asis)
    g2_opt(string asis)
    unstandardized
    graph_export(string asis)
    msize(passthru)
    mattitles
    method(string asis)
    returnweights
    generate(string asis)
    ]
    ;
#delimit cr  


*------------------------------------------------------------------------------*
* (0) Error checks in parsing
*------------------------------------------------------------------------------*
tempvar touse
mark `touse' `if' `in'
**Check if group ID is numeric or string
local clustvar "`2'"

local stringvar=0
cap count if `2'==0
if _rc!=0 {
    local stringvar=1
    local groupvar `2'
    tempvar ID
    egen `ID' = group(`2')
    local varlist `1' `ID' `3' `4'
    tokenize `varlist'
}
else {
    tokenize `varlist'
}

if !inlist("`vce'", "bootstrap", "jackknife", "placebo","noinference") {
    dis as error "vce() must be one of bootstrap, jackknife, placebo or noinference."
    exit 198
}

if (length("`if'")+length("`in'")>0) {
    preserve
    qui keep if `touse'
}

qui xtset `2' `3'
if `"`r(balanced)'"'!="strongly balanced" {
    dis as error "Panel is unbalanced."
    exit 451
}
qui count if `1'==.
if r(N)!=0 {
    dis as error "Missing values found in dependent variable.  A balanced panel without missing observations is required."
    exit 416
}
qui count if `4'==.
if r(N)!=0 {
    dis as error "Missing values found in treatment variable.  A balanced panel without missing observations is required."
    exit 416
}
qui count if `4'!=0&`4'!=1
if r(N)!=0 {
    dis as error "Treatment variable takes values distinct from 0 and 1."
    exit 450
}
qui sum `3'
qui sum `4' if `3'==r(min)
if r(max)!=0 {
    dis as error "Certain units are treated in the first period of the panel."
    dis as error "Any units which are treated the entire panel should not be included in the synthetic DID procedure."
    exit 459
}
qui sum `4'
if (r(min)==0 & r(max)==0)==1 {
    di as error "All units are controls."
    exit 459
}
tempvar a1 a2
qui gen `a1'=`3' if `4'==1
qui bys `2': egen `a2'=min(`a1')
qui levelsof `a2', local(A2)
foreach t of local A2 {
    qui sum `4' if `3'>=`t'
    if r(min)==r(max) {
        di as error "All units are treated in some adoption time."
        exit 459
    }
}

tempvar test
qui bys `2' (`3'): gen `test'=`4'-`4'[_n-1] 
qui sum `test'
qui count if `test'!=0&`test'!=1&`test'!=.
if r(N)!=0 {
    local e1 "to only change from untreated to treated, or remain untreated."
    dis as error "Units are observed to change from treated (earlier) to untreated (later)."
    dis as error "A staggered adoption is assumed in which units are assumed `e1'"
    exit 459
}
drop `test'

local SEN "standard error needs"
if "`vce'"=="jackknife" {
    tempvar t1 t2
    qui gen `t1'=`3' if `4'==1
    qui bys `2': egen `t2'=min(`t1')
    qui levelsof `t2', local(T2)
    foreach t of local T2 {
        qui sum `2' if `4'==1 & `t2'==`t'
        if r(min)==r(max) {
            di as error "Jackknife `SEN' at least two treated units for each treatment period"
            exit 451
        }
    }
}

if "`vce'"=="bootstrap" {
    tempvar t1 t2
    qui gen `t1'=`3' if `4'==1
    qui bys `2': egen `t2'=min(`t1') 
    qui sum `t2' 
    if r(min)==r(max) {
        qui sum `2' if `4'==1
        if r(min)==r(max) {
            local errBS "there is a single treatment period"
            di as error "Bootstrap `SEN' more than one treated unit if `errBS'"
            exit 451
        }
    }
}

if "`vce'"=="placebo" {
    tempvar t1 t2
    qui gen `t1'=`3' if `4'==1 
    qui bys `2': egen `t2'=min(`t1') 
    qui levelsof `t2', local(T2)
	
    local trNtot=0
    foreach t of local T2 {
        qui count if `4'==0 & `3'==`t' & `t2'==.
        local coN=r(N)
        qui count if `4'==1 & `3'==`t' & `t2'==`t'
        local trN=r(N)
		
        local trNtot=`trNtot'+`trN'
        *if (`coN'<`trN')==1 {
         if (`coN'<=`trNtot')==1 {
            di as error "Placebo `SEN' to have more total controls than total treated units"
            exit 451
        }
    }
}

if (length("`if'")+length("`in'")>0) {
    restore
}

*------------------------------------------------------------------------------*
* (1) Basic set-up 
*------------------------------------------------------------------------------*
tempvar treated ty tyear n
qui gen `ty' = `3' if `4'==1 
qui bys `2': egen `treated' = max(`4') if `touse'
qui by  `2': egen `tyear'   = min(`ty') if `touse'

if (length("`if'")+length("`in'")>0) {
    preserve
    qui keep if `touse'
}

sort `3' `treated' `2'
if "`mattitles'"!="" {
    if `stringvar'==1 qui levelsof `groupvar' if `treated'==0
    else              qui levelsof `2' if `treated'==0
    local rnames = r(levels)
}
    
gen `n'=_n
qui sum `3'
local mint=r(min)
qui putmata ori_id=`2' ori_pos=`n' if `3'==`mint' & `tyear'==., replace
if (length("`if'")+length("`in'")>0) {
    restore
}

tempname id_controls
if length("`graph'")!=0&`stringvar'==1 {
    preserve
    qui sum `3' if `touse'
    **Save original state names for later use with graph
    qui keep if `3' == r(min) & `touse'
    keep `groupvar' `2'
    tempfile stateString
    rename `groupvar' stateName
    rename `2' `id_controls'
    qui save "`stateString'"
    restore
}

local control_opt = 0
if "`covariates'"!="" {
    _parse_X `covariates'
    if "`r(ctype)'"=="optimized" local control_opt = 2
    if "`r(ctype)'"=="projected" local control_opt = 1
    local conts = r(controls)
    local contname = r(controls)
	
    foreach var of varlist `contname' {
        qui sum `var'
        if r(sd)==0 {
            dis as error "Covariates were found to be constant in the estimation sample."
            dis as error "Please remove constant covariates from the covariate set."
            exit 416
        }
    }
	
    if `control_opt'==2&length(`"`unstandardized'"')==0 {
        local sconts
        foreach var of varlist `conts' {
            tempvar z`var'
            qui egen `z`var''= std(`var') if `touse'
            local sconts `sconts' `z`var''
        }
        local conts `sconts'
    }
	
    tempvar nmiss
    qui egen `nmiss' = rowmiss(`conts')
    qui sum `nmiss' if `touse'
    if r(mean)!=0 {
        dis as error "Missing values found in covariates."
        dis as error "A balanced panel without missing observations is required."
        exit 416
    }
}

*------------------------------------------------------------------------------*
* (2) Calculate ATT, plus some locals for inference
*------------------------------------------------------------------------------*
if "`vce'"=="jackknife" local jk=1
else local jk=0

if (length("`if'")+length("`in'")>0) {
    preserve
    qui keep if `touse'
}

if "`method'"==""     local m=1
if "`method'"=="sdid" local m=1
if "`method'"=="did"  local m=2
if "`method'"=="sc"   local m=3

**IDs (`2') go in twice here because we are not resampling
mata: data = st_data(.,("`1' `2' `2' `3' `4' `treated' `tyear' `conts'"))
mata: ATT = synthdid(data, 0, ., ., `control_opt', `jk', `m', `zeta_lambda', `zeta_omega', `min_dec', `max_iter')
mata: LAMBDA = st_matrix("lambda")
mata: OMEGA  = st_matrix("omega")

*Warning for reached maximun number of iteration
if (`m'==1) {
    if `reached_lambda'==1 {
        local reached_l "Maximun number of iterations (`max_iter') reached for lambda weight calculation"
    }
}
if (`m'==1 | `m'==3) {
    if `reached_omega'==1 {
        local reached_o "Maximun number of iterations (`max_iter') reached for omega weight calculation"
    }
}
mata: tau    = st_matrix("tau") 
mata: st_local("ATT", strofreal(ATT))

if "`jk'"=="1" {
    mata: se_tau = st_matrix("se_i")
    mata: st_matrix("se_tau",se_tau)
}

qui count if `3'==`mint'                //total units
local N=r(N)
qui count if `treated'==0 & `3'==`mint' //control units
local co=r(N)
qui count if `treated'==1 & `3'==`mint' //treated units (total)
local tr=r(N)
local newtr=`co'-`tr'+1                 //start of treated units
qui levelsof `3', local(ttime)          //T
local T: word count `r(levels)'
local Ntotal=`N'*`T'                    //total observations

matrix ttime = J(`T',1,.)
local jj=1
foreach l of local ttime {
    matrix ttime[`jj',1]=`l'
    local ++jj
}

mkmat `tyear' if `tyear'!=. & `3'==`mint', matrix(tryears) //save adoption values

qui levelsof `tyear', local(tryear) //adop
local nadop: word count `r(levels)'
matrix adoption = J(`nadop',1,.)
local jj=1
foreach l of local tryear {
    matrix adoption[`jj',1]=`l'
    local ++jj
}

if (length("`if'")+length("`in'")>0) {
    restore
}

if "`vce'"!="noinference" {
*--------------------------------------------------------------------------*
* (3) Standard error: bootstrap
*--------------------------------------------------------------------------*
if "`vce'"=="bootstrap" {
    cap set seed `seed'
    local b = 1
    local B = `reps'
    mata: ATT_b = J(`B', 1, .)
    mata: tau_b = J(1, 2, .)
    mata: se_tau = J(`nadop', 1, .)
	
    dis "Bootstrap replications (`reps'). This may take some time."
    dis "----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5"
	
    while `b'<=`B' {
        preserve
        qui keep if `touse'
        keep `1' `2' `3' `4' `treated' `tyear' `conts'
        tempvar cID
        bsample, cluster(`2') idcluster(`cID')
        qui count if `treated' == 0
        local r1 = r(N)
        qui count if `treated' != 0
        local r2 = r(N)
        
        if (`r1'==0 | `r2'==0) {
            *di "all units are control or treated"
        }	
        else {
            display in smcl "." _continue
            if mod(`b',50)==0 dis "     `b'"

            qui putmata bsam_id=`2' if `tyear'==. & `3'==`mint', replace
            mata: indicator=smerge(bsam_id, (ori_id, ori_pos))
            mata: OMEGAB=OMEGA[(indicator\rows(OMEGA)),]		
            mata: data = st_data(.,("`1' `cID' `2' `3' `4' `treated' `tyear' `conts'"))
            mata: ATT_b[`b',] = synthdid(data,1,OMEGAB,LAMBDA,`control_opt',`jk',`m', `zeta_lambda', `zeta_omega', `min_dec', `max_iter')
			
            *taus for adoption times
            mata: ty = uniqrows(select(data[,7], data[,5]:==1))
            mata: ty_taus = (ty , st_matrix("tau_i"))
            mata: tau_b =  (tau_b\ty_taus)
			
            local ++b
        }
        restore
    }
	
    *SE for adoption times
    local i=1
    foreach tr of local tryear {
        mata: nr_tau = rows(select(tau_b[,2], tau_b[,1]:==`tr'))
        mata: vect_tau = select(tau_b[,2], tau_b[,1]:==`tr')
        mata: se_tau[`i',1] = sqrt((nr_tau-1)/nr_tau) * sqrt(variance(vec(vect_tau)))
        mata: st_matrix("se_tau",se_tau)	
        local ++i
    }
	
    mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(ATT_b)))
    mata: st_local("se", strofreal(se)) 
}
*--------------------------------------------------------------------------*
* (4) Standard error: placebo
*--------------------------------------------------------------------------*
else if "`vce'"=="placebo" {
    cap set seed `seed'
    local b = 1
    local B = `reps'
    mata: ATT_p = J(`B', 1, .)
    mata: tau_p = J(1, 2, .)
    mata: se_tau = J(`nadop', 1, .)
	
    dis "Placebo replications (`reps'). This may take some time."
    dis "----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5"
	
    while `b'<=`B' {
        preserve
        qui keep if `touse'
		
        keep `1' `2' `3' `4' `tyear' `conts'
        qui drop if `tyear'!=.        //drop treated units
		
        tempvar id
        egen `id' = group(`2')
        qui sum `id'
        local max = r(max)
		
        *generate vector of selected placebos
        local rowsselect = `co'-`newtr'+1				
        mata id_list = range(1, `max', 1)
        mata asignar_pr = (id_list, runiform(`max',1))
        mata st_matrix("Select", sort(asignar_pr, 2)[1..`rowsselect', 1])
        				
        *replace treatment status with placebos
        forval i=1/`rowsselect' {
            qui replace `tyear' = tryears[`i',1] if `id'==Select[`i',1]
        }

        qui replace `4' = 1 if `3'>=`tyear'
        egen `treated' = max(`4'), by(`2')
				
        display in smcl "." _continue
        if mod(`b',50)==0 dis "     `b'"
		
        qui putmata psam_id=`2' if `tyear'==. & `3'==`mint', replace
        mata: indicator=smerge(psam_id, (ori_id, ori_pos))
        mata: OMEGAP=OMEGA[(indicator\rows(OMEGA)),]
        mata: data = st_data(.,("`1' `2' `2' `3' `4' `treated' `tyear' `conts'"))
        mata: ATT_p[`b',] = synthdid(data,1,OMEGAP,LAMBDA,`control_opt',`jk',`m', `zeta_lambda', `zeta_omega', `min_dec', `max_iter')

        *taus for adoption times
        mata: ty = uniqrows(select(data[,7], data[,5]:==1))
        mata: ty_taus = (ty, st_matrix("tau_i"))
        mata: tau_p =  (tau_p\ty_taus)
							
        local ++b
        restore
    }
	
    *SE for adoption times
    local i=1
    foreach tr of local tryear {
        mata: nr_tau = rows(select(tau_p[,2], tau_p[,1]:==`tr'))
        mata: vect_tau = select(tau_p[,2], tau_p[,1]:==`tr')
        mata: se_tau[`i',1] = sqrt((nr_tau-1)/nr_tau) * sqrt(variance(vec(vect_tau)))
        mata: st_matrix("se_tau",se_tau)	
        local ++i
    }
    mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(ATT_p)))
    mata: st_local("se", strofreal(se))
}
}
*--------------------------------------------------------------------------*
* (5) Return output
*--------------------------------------------------------------------------*
ereturn clear

*for esttab and eststo
if "`vce'"!="noinference" {
matrix b=`ATT'
matrix V=`se'^2
matrix colnames b=`4'
matrix rownames b=`1'
matrix colnames V=`4'
matrix rownames V=`4'
ereturn post b V, depname(`1') obs(`Ntotal')
}

if "`vce'"!="noinference" {
	matrix tau=(tau,se_tau,adoption)
	matrix colnames tau = Tau Std.Err. Time
}

else {
	matrix tau=(tau,adoption)
	matrix colnames tau = Tau Time
}

qui levelsof `2' if `touse'
local nclust: word count `r(levels)'

if "`vce'"=="noinference" local se=.
ereturn scalar se     =`se' 
ereturn scalar ATT    =`ATT'
ereturn scalar reps   =`reps'
ereturn scalar N_clust=`nclust'

if "`mattitles'"!="" {
        local rn ""
        foreach n of local rnames {
            local n = subinstr("`n'", "`", " ", .)
            local n = subinstr("`n'", "'", " ", .)
            local rn `rn' `=strtoname("`n'")'
        }
        local rn `rn' "Adoption Time"
        matrix rownames omega = `rn'
}

if length("`returnweights'")!=0 {
    preserve
    clear 
    qui svmat lambda, n(col)
	
    local i=1
    foreach t of local tryear {
	    ren c`i' `generate'lambda`t'
        if (`m'==1 | `m'==2) {
            qui count if `generate'lambda`t'==.
            qui replace `generate'lambda`t' = 1/r(N) if `generate'lambda`t'==.
        }
        else {
            qui count if `generate'lambda`t'==0
            qui replace `generate'lambda`t' = 1/r(N) if `generate'lambda`t'==0
            qui count if `generate'lambda`t'==.
            qui replace `generate'lambda`t' = 1/r(N) if `generate'lambda`t'==.
        }
        local ++i
    }
    ren c `3'
    tempfile dlambda
    qui save `dlambda'
    restore

    qui merge m:1 `3' using `dlambda', nogen keep(1 3)

    preserve
    clear 
    qui svmat omega, n(col)
	
    local i=1
    foreach t of local tryear {
	    ren c`i' `generate'omega`t'
        local ++i
    }
    ren c `2'
    tempfile domega
    qui save `domega'
    restore

    qui merge m:1 `2' using `domega', nogen keep(1 3)
	
    local i=1
    foreach t of local tryear {
        qui replace `generate'omega`t' = 1/trunitN[`i',1] if `generate'omega`t'==.
        local ++i
    }
	
}

ereturn matrix tau      tau
ereturn matrix lambda   lambda
ereturn matrix omega    omega
ereturn matrix adoption adoption


if "`covariates'"!="" {
    if "`control_opt'"=="1" matrix rownames beta = `contname'
    if "`control_opt'"=="2" matrix rownames beta = `contname' adoption
    ereturn matrix beta beta
}

ereturn local cmdline  "sdid `0'"
ereturn local clustvar `clustvar'
ereturn local depvar   `1'
ereturn local vce      "`vce'"
ereturn local title    "Synthetic Difference-in-Differences"
ereturn local cmd      "sdid"

local t=`ATT'/`se'
local pval= 2 * (1-normal(abs(`ATT'/`se')))
local v = (1-`level'/100)/2
local LCI = `ATT' + invnormal(`v')*`se'
local UCI = `ATT' + invnormal(1-`v')*`se'
ereturn scalar ATT_l = `LCI' 
ereturn scalar ATT_r = `UCI'

if ("`method'"=="sdid" | "`method'"=="" ) {
    local tabletitle "Synthetic Difference-in-Differences Estimator"
    local tablefootnote "Refer to Arkhangelsky et al., (2020) for theoretical derivations."
}
if ("`method'"=="did") {
    local tabletitle "Difference-in-Differences Estimator"
    local tablefootnote
}
if ("`method'"=="sc") {
    local tabletitle "Synthetic Control"
    local tablefootnote
}


di as text ""
di as text ""
di as text "`tabletitle'" 
di as text ""
di as text "{hline 13}{c TT}{hline 63}"
di as text %12s abbrev("`1'",12) " {c |}     ATT     Std. Err.     t      P>|t|    [`level'% Conf. Interval]" 
di as text "{hline 13}{c +}{hline 63}"
di as text %12s abbrev("`4'",12) " {c |} " as result %9.5f `ATT' "  " %9.5f `se' %9.2f `t' %9.3f `pval' "   " %9.5f `LCI' "   " %9.5f `UCI'
di as text "{hline 13}{c BT}{hline 63}"
di as text "95% CIs and p-values are based on Large-Sample approximations."
di as text "`tablefootnote'" 
di as text "`reached_l'"
di as text "`reached_o'"

*--------------------------------------------------------------------------*
* (6) Graphing
*--------------------------------------------------------------------------*
    tempname time lambdadrop
    matrix auxlambda = e(lambda)
    qui svmat series
    ren series1 `time'
    la var `time' `3'
    qui svmat auxlambda
	
    local ColConTr 
    local ColConTrDiff 
    local colcon = 2
    local coltr  = 3
    local colaux = 1
    local colauxfin = 1
    foreach t of local tryear {
        tempname Yco`t' Ytr`t' lambda`t'
        local ColConTr `ColConTr' Yco`t' Ytr`t'
        local ColConTrDiff `ColConTrDiff' Diff`t'
        ren series`colcon'    `Yco`t''
        ren series`coltr'     `Ytr`t''
        ren auxlambda`colaux' `lambda`t''

        local colcon    = `colcon' + 2
        local coltr     = `coltr'  + 2
        local colaux    = `colaux' + 1
        local colauxfin = `colauxfin' + 1
    }	
    ren auxlambda`colauxfin' `lambdadrop'	
	
    mat coln series = `3' `ColConTr'
    mat coln difference = `3' `ColConTrDiff'
    ereturn matrix series series
    ereturn matrix difference difference
	
    if length("`graph'")!=0 {
        if (`co'>1000 & "`g1on'"=="g1on") {
            local vis "and graphed unit weights are unlikely to be easily visualized."
            dis "Be careful with graph option. You have a lot of control units `vis'"
        }

        if length(`"`graph_export'"')!=0 {
            _graph_Name `graph_export'
            local gstub = r(gstub)
            local suffix = r(suffix)
            			
            if "`suffix'"==".gph" local pre graph save
            else                  local pre graph export
			
            if "`gstub'"=="."  local gstub
            if "`suffix'"=="." local suffix
            local ex=1
        }
        else local ex=0

        if "`g1on'"=="g1on" {
            tempvar order
            matrix aux_tau = e(tau)
            qui svmat dif_plot
            ren dif_plot1 `id_controls'
            qui gen `order' = _n if `id_controls'!=.
			
            if `stringvar'== 0 {
                qui levelsof `id_controls', local(sgroup)
                foreach s of local sgroup {
                    qui sum `order' if `id_controls'==`s'
                    local oN = r(mean)
                    if r(mean)!= . {
                        local xlabs `xlabs' `oN' "`s'"
                    }
                }
            }

            if `stringvar'==1 {
                qui merge m:1 `id_controls' using "`stateString'", keep(1 3) nogen
                qui levelsof stateName, local(sgroup)

                foreach s of local sgroup {
                    qui sum `order' if stateName == `"`s'"'
                    local oN = r(mean)
                    if r(mean)!= . {
                        local xlabs `xlabs' `oN' "`s'"
                    }
                }            
            }

            local ColDif1 = 2
            local ColDif2 = 3
            local aux_tau = 1
            foreach t of local tryear {
                tempname Diff`t' Weight`t'
                ren dif_plot`ColDif1' `Diff`t''
                ren dif_plot`ColDif2' `Weight`t''

                if length(`"`msize'"')==0 local ms msize(tiny)
                else local ms `msize'
                lab var `Diff`t'' "Difference"
                lab var `order' "Group"		
                local tau = aux_tau[`aux_tau',1]
				
                #delimit ;
                twoway scatter `Diff`t'' `order' if `Weight`t''!=0 & `id_controls'!=. [aw=`Weight`t''], `ms'
                    || scatter `Diff`t'' `order' if `Weight`t''==0 & `id_controls'!=., m(X) 
                    xlabel(`xlabs', angle(vertical) labs(vsmall) valuelabel)
                    yline(`tau', lc(red))
                    `g1_opt' name(g1_`t', replace) legend(off);
                if `ex'==1 `pre' "`gstub'weights`t'`suffix'", replace;
                #delimit cr
				
                local ColDif1 = `ColDif1' + 2
                local ColDif2 = `ColDif2' + 2
                local ++aux_tau
            }
            if `stringvar'==1 cap drop stateName
        }

        foreach t of local tryear {
            local tit   "ytitle("Lambda weight", margin(0 0 0 45) axis(2))"
            local ylab  "ylabel(0(0.25)1, format(%5.1f) axis(2))"
            local scale "yscale(range(0(1)3) axis(2))"
            local axis  "yaxis(2) lp(solid)"
            local timelambda "|| area `lambda`t'' `time' if `time'!=., `axis' `ylab' `scale' `tit'"
            if "`method'"=="sc" local timelambda 
		
            #delimit ;
            twoway line `Yco`t'' `time', yaxis(1) lp(dash)
                || line `Ytr`t'' `time', yaxis(1) lp(solid)
                `timelambda' ||, 
            xline(`t', lc(red)) legend(order(1 "Control" 2 "Treated") pos(12) col(2))
            `g2_opt' name(g2_`t', replace);
            if `ex'==1 `pre' "`gstub'trends`t'`suffix'", replace;
            #delimit cr
        }
    }
	
end

*------------------------------------------------------------------------------*
* (7) Stata Subroutines
*------------------------------------------------------------------------------*
cap program drop _parse_X
program define _parse_X, rclass
    syntax varlist [, *]
    
    local valid_X=inlist("`options'", "optimized", "projected", "")
    if (`valid_X'==0) {
        dis as error "Control types must be one of {bf:optimized} or {bf:projected}."
        exit 198
    }
    
    if ("`options'"=="") local options optimized
    unab conts: `varlist'
    
    return local controls `conts'
    return local ctype    `options'
end

cap program drop _graph_Name
program define _graph_Name, rclass
    syntax [anything] [, *]
    
    local valid_X=inlist("`options'", ".ps", ".eps", ".svg", ".wmf", ".emf", ".pdf", ".png", ".tif", ".gph")
    if (`valid_X'==0) {
        dis as error "When specifying graph export, export format must be a valid type."
        dis as error "Please specify a valid export type (.eps, .ps, .pdf, and so forth)"
        exit 198
    }

    return local gstub  `anything'
    return local suffix `options'
end

*------------------------------------------------------------------------------*
* (8) Mata functions
*------------------------------------------------------------------------------*
*Main function (synthdid)
*Below assumes data is a matrix with:
* (1) y variable
* (2) group variable (bootstrap fix if necessary)
* (3) orig group variable
* (4) time variable
* (5) treatment variable
* (6) indicator if unit ever treated
* (7) indicator of year treated (if treated)
* (8+) any controls
mata:
    real scalar synthdid(matrix data, inference, OMEGA, LAMBDA, controls, jk, mt, ELambda, EOmega, MinDec, MaxIter) {
        data  = sort(data, (6,2,4))
        units = panelsetup(data,2)
        NT = panelstats(units)[,3]
        treat=panelsum(data[.,(2,5)],units)
        treat[,1]=treat[,1]/NT
        treat[,2]=1*(treat[,2]:>=1) + 0*(treat[,2]:==1)
        Nco = sum(data[,7]:==.)/NT
        controlID = uniqrows(select(data[.,2],data[,7]:==.))  
		
        if (jk==1) {
            uniqID=(uniqrows(select(data[.,2],data[,7]:==.)) \ uniqrows(select(data[.,2],data[,7]:!=.)))
            N = panelstats(units)[,1]
        }
		
        //matrix for beta
        if (inference==0 & (controls==0)) {
            Beta = J(1, 1, .)
        }
		
        //matrix for control, treated and difference series (graph option)
        Series = J(NT,1,.)
        Diff = J(NT,1,.)
        D_plot = J(Nco,1,.)
		
        //Adjust for controls in xysnth way
        //save original data for jackknife
        data_ori=data
        if (cols(data)>7 & controls==1) {
            projected(data, Yprojected, Beta)
            data[,1] = Yprojected
        }
        
        //Find years which change treatment
        trt = select(uniqrows(data[,7]),uniqrows(data[,7]):!=.)
		
        //Iterate over years, calculating each estimate
        tau    = J(rows(trt),1,.)
        tau_wt = J(1,rows(trt),.)
		
        //Treated unit vector
        trunitN = J(rows(trt),1,.)
		
        if (inference==0) {
            Omega = J(Nco, rows(trt),.)
            Lambda = J(NT, rows(trt),.)
        }
        if (inference==0 & controls==2) {
            Beta = J(cols(data)-7, rows(trt),.)
        }
        for(yr=1;yr<=rows(trt);++yr) {
            cond1 = data[,7]:==trt[yr]
            cond2 = data[,7]:==.
            cond = cond1+cond2
            yNtr = sum(cond1)
            yNco = sum(cond2)
            ydata = select(data,cond)
            yunits = panelsetup(ydata,2)
            yNT = panelstats(yunits)[,3]
            yNG = panelstats(yunits)[,1]
            yN  = panelstats(yunits)[,2]
            yNtr = yNtr/yNT
            trunitN[yr] = yNtr 
            yNco = yNco/yNT
            yTpost = max(ydata[,4])-trt[yr]+1
            pretreat = select(uniqrows(data[,4]),uniqrows(data[,4]):<trt[yr])
            Npre  = rows(pretreat)
            Npost = yNT-Npre
            //Calculate Zeta
            ndiff = yNG*(yNT-1)
            ylag = ydata[1..(yN-1),1]
            ylag = (. \ ylag)
            diff = ydata[.,1]-ylag
            first = mod(0::(yN-1),yNT):==0
            postt = ydata[,4]:>=trt[yr]
            dropc = first+postt+ydata[,6]
            prediff = select(diff,dropc:==0)
            sig_t = sqrt(variance(prediff))
            EtaLambda = ELambda
			
            if (mt!=3) {
                if (EOmega==1e-6) {
                    EtaOmega = (yNtr*yTpost)^(1/4)
                }
                else {
                    EtaOmega = EOmega
                }
            }
            if (mt==3) {
                EtaOmega = EOmega
            }
			
            yZetaOmega  = EtaOmega*sig_t
            yZetaLambda = EtaLambda*sig_t
            //Generate Y matrices
            ids = uniqrows(ydata[.,2])
            //ids_trt = uniqrows(select(ydata[.,2], ydata[,7]:==trt[yr]))
            ids_ctr = uniqrows(select(ydata[.,2], ydata[,7]:==.))
            ytreat = ydata[,7]:==trt[yr]
            ytreat=panelsum(ytreat,yunits)
            ytreat=ytreat/yNT
            Y = rowshape(ydata[.,1],yNG)
            Y0 = select(Y,ytreat:==0)
            Y1 = select(Y,ytreat:==1)
            //Generate average of outcomes for each unit over time
            promt = mean(Y0[,(Npre+1)::yNT]')'
            //Calculate input matrices (pre-treatment and averages)
            A_l = Y0[,1..Npre]:-mean(Y0[,1..Npre])
            b_l = promt:-mean(promt)
			
            if (mt!=3) {
                A_o = Y0[,1..Npre]':-mean(Y0[,1..Npre]')
                b_o = mean(Y1[.,1..Npre])':-mean(mean(Y1[.,1..Npre])')
            }
            if (mt==3) {
                A_o = Y0[,1..Npre]'
                b_o = mean(Y1[.,1..Npre])'
            }			
	
            //Calculate Tau for t
            if (inference==1 & mt!=2) {
                lambda_l = select(LAMBDA'[.,1::Npre],LAMBDA'[,rows(LAMBDA)]:==trt[yr])
                l_o = select(OMEGA'[.,1::yNco],OMEGA'[,rows(OMEGA)]:==trt[yr])
                lambda_o = sum_norm(l_o) //update so prior weights sum to 1
            }
			
            else {
                if (mt==3) {
                    lambda_l = J(1,cols(A_l),0)
                }
                if (mt!=3) {
                    lambda_l = J(1,cols(A_l),1/cols(A_l))
                }
                lambda_o = J(1,cols(A_o),1/cols(A_o))
            }
            mindec = (MinDec*sig_t)^2
			
            if (controls==2) {
                K = cols(data)
                A = J((K-7),1,NULL)
                CX = J((K-7),1,NULL)
                n=1
                for (c=8;c<=K;++c) {
                    X = rowshape(ydata[.,c],yNG)
                    CX[n] = &(rowshape(ydata[.,c],yNG))
                    X0 = select(X,ytreat:==0)
                    X1 = select(X,ytreat:==1)
                    A[n] = &((X0[,1..Npre],mean(X0[,(Npre+1)::yNT]')')\(mean(X1[.,1..Npre]),0))
                    n++
                }
				
                maxiter = MaxIter
                vals = J(1, maxiter, .)
                beta=J(1,(K-7),0)
                gradbeta = J(1,(K-7),0)
                t=0
                dd=1
				
                eta_o = Npre*yZetaOmega^2
                eta_l = yNco*yZetaLambda^2
			
                //update wights: sdid and sc
                if (mt!=2) {
                    if (mt==1) {
                        lambda_l = fw(A_l,b_l,lambda_l,eta_l)
                    }
                    lambda_o = fw(A_o,b_o,lambda_o,eta_o)
                }
                err_l    = (A_l, b_l)*(lambda_l' \ -1)
                err_o    = (A_o, b_o)*(lambda_o' \ -1)
				
                while (t<maxiter & (t<2 | dd>mindec)) {
                    t++		
                    for (c=1;c<=(K-7);++c) {
                        gradbeta[c]=-(err_l'*((*A[c])[1..yNco,1..(Npre+1)])*((lambda_l'\-1):/yNco) +
                            err_o'*((*A[c])[1..(yNco+1),1..Npre])'*((lambda_o'\-1):/Npre))
                    }			
                    alpha=1/t
                    beta=beta-alpha*gradbeta
                    //~ contract3
                    C = J(yNco+1,Npre+1,0)
                    for (c=1;c<=(K-7);++c) {
                        C = C + beta[c]*(*A[c])
                    }
                    Ybeta=((Y0[,1..Npre],promt)\(mean(Y1[.,1..Npre]),0))-C
					
                    Ybeta_A_l = Ybeta[1::yNco,1::Npre]:-mean(Ybeta[1::yNco,1::Npre])
                    Ybeta_b_l = Ybeta[1::yNco,Npre+1]:-mean(Ybeta[1::yNco,Npre+1])
                    if (mt!=3) {
                       Ybeta_A_o = Ybeta[1::yNco,1::Npre]':-mean(Ybeta[1::yNco,1::Npre]')
                       Ybeta_b_o = (Ybeta[yNco+1,1::Npre])':-mean((Ybeta[yNco+1,1::Npre])')
                    }
                    if (mt==3) {
                       Ybeta_A_o = Ybeta[1::yNco,1::Npre]'
                       Ybeta_b_o = (Ybeta[yNco+1,1::Npre])'
                    }
				
                    if (mt!=2) {
                        if (mt==1) {
                            lambda_l = fw(Ybeta_A_l,Ybeta_b_l,lambda_l,eta_l)
                        }
                        lambda_o = fw(Ybeta_A_o,Ybeta_b_o,lambda_o,eta_o)
                    }
                    err_l    = (Ybeta_A_l,Ybeta_b_l)*(lambda_l' \ -1)
                    err_o    = (Ybeta_A_o,Ybeta_b_o)*(lambda_o' \ -1)
                    vals[1,t]=yZetaOmega^2*(lambda_o*lambda_o')+yZetaLambda^2*(lambda_l*lambda_l')+
                              (err_o'*err_o)/Npre+(err_l'*err_l)/yNco
                    if (t>1) dd = vals[1,t-1] - vals[1,t]
                }
                if (inference==0) {
                    Lambda[.,yr] =  (lambda_l' \ J(Npost,1,.))
                    Omega[.,yr] = lambda_o'
                    Beta[.,yr] = beta'
                }
                Xbeta = J(rows(Y),cols(Y),0)
                for (c=1;c<=(K-7);++c) {
                    Xbeta = Xbeta + beta[c]*(*CX[c])
                }
                Y=Y-Xbeta
                tau[yr] = (-lambda_o, J(1,yNtr,1/yNtr))*(Y)*(-lambda_l, J(1,Npost,1/Npost))'
                tau_wt[yr] = yNtr*Npost
				
                //for line graph
                Diff = (Diff, ((-lambda_o, J(1,yNtr,1/yNtr))*Y)')
                Series = (Series, ((lambda_o)*Y0)', ((J(1,yNtr,1/yNtr))*Y1)')
                //for scatter graph	
                Npre_min = Npre + 1
                D = (mean(Y[,Npre_min::cols(Y)]')' - (lambda_l*Y[,1::Npre]')', (J(yNco,1,0)\J(yNtr,1,1)))
                Daux = mean(select(D[,1], D[,2]:==1))
                D_plot = (D_plot, ((Daux:-D[,1]), (lambda_o'\J(yNtr,1,.)))[1::yNco,])
            }
			
            if (controls==0 | controls==1) {
                if (mt!=2) {
                    //Find optimal weight matrices
                    eta_o = Npre*yZetaOmega^2
                    eta_l = yNco*yZetaLambda^2
                    if (mt==1) {
                        lambda_l = lambda(A_l,b_l,lambda_l,eta_l,yZetaLambda,100,mindec,0)
                        lambda_l = sspar(lambda_l)
                        lambda_l = lambda(A_l, b_l, lambda_l,eta_l,yZetaLambda, MaxIter,mindec,1)
                        st_local("reached_lambda", st_local("reached"))
                    }
                    lambda_o = lambda(A_o, b_o, lambda_o,eta_o,yZetaOmega, 100,mindec,0)
                    lambda_o = sspar(lambda_o)
                    lambda_o = lambda(A_o, b_o, lambda_o,eta_o,yZetaOmega, MaxIter,mindec,1)
                    st_local("reached_omega", st_local("reached"))
                }
                if (inference==0) {
                    Lambda[.,yr] =  (lambda_l' \ J(Npost,1,.))
                    Omega[.,yr] = lambda_o'
                }
                tau[yr] = (-lambda_o, J(1,yNtr,1/yNtr))*Y*(-lambda_l, J(1,Npost,1/Npost))'
                tau_wt[yr] = yNtr*Npost
				
                //mu = (-lambda_o, J(1,yNtr,1/yNtr))*Y*(lambda_l, J(1,Npost,0))'

                //for line graph
                Diff = (Diff, ((-lambda_o, J(1,yNtr,1/yNtr))*Y)')
                Series = (Series, ((lambda_o)*Y0)', ((J(1,yNtr,1/yNtr))*Y1)')
                //for scatter graph	
                Npre_min = Npre + 1
                D = (mean(Y[,Npre_min::cols(Y)]')' - (lambda_l*Y[,1::Npre]')', (J(yNco,1,0)\J(yNtr,1,1)))
                Daux = mean(select(D[,1], D[,2]:==1))
                D_plot = (D_plot, ((Daux:-D[,1]), (lambda_o'\J(yNtr,1,.)))[1::yNco,])

            }
        }
		
        if (inference==0) {
            Omega =  (Omega, controlID)
            Omega =  (Omega \ (trt', .))
            Lambda = (Lambda, uniqrows(data[,4]))
            Lambda = (Lambda \ (trt', .))
            Series[,1] = uniqrows(data[,4])
            Diff[,1] = uniqrows(data[,4])
            D_plot[,1] = (ids_ctr) 
        }
        if (inference==0 & controls==2) {
            Beta = (Beta \ trt')
        }		
		
        //jackknife
        if (jk==1) {
            tau_aux    = J(rows(trt),1,.)
            tau_wt_aux = J(1,rows(trt),.)
            ATT_aux = J(N,1,.)
            tau_aux_j = J(rows(trt),N,.)
			se_aux_j = J(rows(trt),1,.)
            yNco_ori=yNco
            ind=(1::yNco)
            for (i=1; i<=N; ++i) {
                drp=uniqID[i]
                data_aux = select(data_ori, data_ori[,3]:!=drp)
                //projected adjustment
                if (cols(data_aux)>7 & controls==1) {
                    projected(data_aux, Yprojected, Beta_jk)
                    data_aux[,1] = Yprojected
                }	
                
                for (yr=1;yr<=rows(trt);++yr) {
                    cond1 = data_aux[,7]:==trt[yr]
                    cond2 = data_aux[,7]:==.
                    cond = cond1+cond2
                    yNtr = sum(cond1)
                    yNco = sum(cond2)
                    yNco = yNco/yNT
                    yNtr = yNtr/yNT
                    pretreat = select(uniqrows(data_aux[,4]),uniqrows(data_aux[,4]):<trt[yr])
                    Npre  = rows(pretreat)		
                    Npost = yNT-Npre
                    ydata_aux = select(data_aux,cond)
                    yunits = panelsetup(ydata_aux,2)
                    yNG = panelstats(yunits)[,1]
                    ytreat = ydata_aux[,7]:==trt[yr]
                    ytreat = panelsum(ytreat,yunits)
                    ytreat = ytreat/yNT
                    Y_aux = rowshape(ydata_aux[.,1],yNG)
					
                    lambda_aux = select(Lambda'[.,1::Npre],Lambda'[,rows(Lambda)]:==trt[yr])  
                    id1=select(ind, ind:!=i)
                    omega_aux = select(Omega'[.,1::yNco_ori],Omega'[,rows(Omega)]:==trt[yr]) 
                    omega_aux = sum_norm(omega_aux[id1])
                    if (controls==0 | controls==1) {
                        tau_aux[yr] = (-omega_aux, J(1,yNtr,1/yNtr))*Y_aux*(-lambda_aux,
                            J(1,Npost,1/Npost))'
                        tau_wt_aux[yr] = yNtr*Npost	
                    }
                    //R adjustment
                    if (cols(data_aux)>7 & controls==2) {
                        maxiter=10000
                        Y0_aux = select(Y_aux,ytreat:==0)
                        Y1_aux = select(Y_aux,ytreat:==1)
                        promt_aux = mean(Y0_aux[,(Npre+1)::yNT]')'
                        A_l_aux = Y0_aux[,1..Npre]:-mean(Y0_aux[,1..Npre])
                        b_l_aux = promt_aux:-mean(promt_aux)
						
                        if (mt!=3) {
                           A_o_aux = Y0_aux[,1..Npre]':-mean(Y0_aux[,1..Npre]')
                           b_o_aux = mean(Y1_aux[.,1..Npre])':-mean(mean(Y1_aux[.,1..Npre])')
                        }
                        if (mt==3) {
                           A_o_aux = Y0_aux[,1..Npre]'
                           b_o_aux = mean(Y1_aux[.,1..Npre])'
                        }
					
                        K = cols(data_aux)
                        A = J((K-7),1,NULL)
                        CX = J((K-7),1,NULL)
                        n=1
                        for (c=8;c<=K;++c) {
                            X = rowshape(ydata_aux[.,c],yNG)
                            CX[n] = &(rowshape(ydata_aux[.,c],yNG))
                            X0 = select(X,ytreat:==0)
                            X1 = select(X,ytreat:==1)
                            A[n] = &((X0[,1..Npre],mean(X0[,(Npre+1)::yNT]')')\(mean(X1[.,1..Npre]),0))
                            n++
                        }
				
                        vals = J(1, maxiter, .)
                        beta=J(1,(K-7),0)
                        gradbeta = J(1,(K-7),0)
                        t=0
                        dd=1	
						
                        //update weights
                        err_l = (A_l_aux, b_l_aux)*(lambda_aux' \ -1)
                        err_o = (A_o_aux, b_o_aux)*(omega_aux' \ -1)	
						
                        while (t<maxiter & (t<2 | dd>mindec)) {
                            t++	
                            for (c=1;c<=(K-7);++c) {
                                gradbeta[c]=
                                    -(err_l'*((*A[c])[1..yNco,1..(Npre+1)])*((lambda_aux'\-1):/yNco) +
                                    err_o'*((*A[c])[1..(yNco+1),1..Npre])'*((omega_aux'\-1):/Npre))
                            }							
                            alpha=1/t
                            beta=beta-alpha*gradbeta
                            C = J(yNco+1,Npre+1,0)
                            for (c=1;c<=(K-7);++c) {
                                C = C + beta[c]*(*A[c])
                            }
                            Ybeta=((Y0_aux[,1..Npre],promt_aux)\(mean(Y1_aux[.,1..Npre]),0))-C
                            Ybeta_A_l = Ybeta[1::yNco,1::Npre]:-mean(Ybeta[1::yNco,1::Npre])
                            Ybeta_b_l = Ybeta[1::yNco,Npre+1]:-mean(Ybeta[1::yNco,Npre+1])
                            Ybeta_A_o = Ybeta[1::yNco,1::Npre]':-mean(Ybeta[1::yNco,1::Npre]')
                            Ybeta_b_o = (Ybeta[yNco+1,1::Npre])':-mean((Ybeta[yNco+1,1::Npre])')
                            err_l = (Ybeta_A_l,Ybeta_b_l)*(lambda_aux' \ -1)
                            err_o = (Ybeta_A_o,Ybeta_b_o)*(omega_aux' \ -1)
							
                            vals[1,t]=yZetaOmega^2*(lambda_aux*lambda_aux')+
                                     yZetaLambda^2*(lambda_aux*lambda_aux')+
                                     (err_o'*err_o)/Npre+(err_l'*err_l)/yNco
                            if (t>1) dd = vals[1,t-1] - vals[1,t]
                        }
						
                        Xbeta = J(rows(Y_aux),cols(Y_aux),0)
                        for (c=1;c<=(K-7);++c) {
                            Xbeta = Xbeta + beta[c]*(*CX[c])
                        }
                        Y_aux=Y_aux-Xbeta
						
                        tau_aux[yr] = (-omega_aux, J(1,yNtr,1/yNtr))*(Y_aux)*(-lambda_aux, J(1,Npost,1/Npost))'
                        tau_wt_aux[yr] = yNtr*Npost	
                    }				
                }
                tau_aux_j[,i] = tau_aux
                tau_wt_aux = tau_wt_aux/sum(tau_wt_aux')
                ATT_aux[i] = tau_wt_aux*tau_aux	
            }
            // SE for any adoption
            for (i=1; i<=rows(trt); ++i) {
                se_aux_j[i,] = sqrt(((N-1)/N) * (N - 1) * variance(vec(tau_aux_j[i,])))
            }
			
            se_j = sqrt(((N-1)/N) * (N - 1) * variance(vec(ATT_aux)))
            st_local("se", strofreal(se_j))
            st_matrix("se_i", se_aux_j)
        }
		
        tau_wt = tau_wt/sum(tau_wt')
        ATT = tau_wt*tau
		
        if (inference==0) {
            st_matrix("omega" ,Omega)
            st_matrix("lambda",Lambda)
            st_matrix("tau"   ,tau)
            st_matrix("beta"  ,Beta)	
            st_matrix("series",Series)
            st_matrix("difference",Diff)
            st_matrix("dif_plot",D_plot)
            st_matrix("trunitN",trunitN)
        }
        if (inference==1 & jk!=1) {
            tau_i = tau
            st_matrix("tau_i",tau_i)
        }
        return(ATT)
    }
end
  
*minimization
mata:
real vector lambda(matrix A, matrix b, matrix x, eta, zeta, maxIter, mindecrease, rev) {    
    row = rows(A)
    col = cols(A)
    vals = J(1, maxIter, .)
    t=0
    dd=1
    
    while (t<maxIter & (t<2 | dd>mindecrease)) {
        t++
        if (rev==1) {
            if (t==maxIter) {
                st_local("reached", strofreal(1))
            }
            else {
                st_local("reached", strofreal(0))
            }
        }
		
        Ax = A * x'	
        hg = (Ax - b)' * A + eta * x
        i = select((1..cols(hg)), colmin(hg :== min(hg)))[1,1]
        dx = -x
        dx[1,i] = 1 - x[1,i]
        v = abs(min(dx))+abs(max(dx))
        if (v==0) {
            x = x
            err = (A, b) * (x' \ -1)
            vals[1,t] = zeta^2 * (x * x') + (err' * err) / row
            if (t>1) {
                dd = vals[1,t-1] - vals[1,t]
            }
        }
        else {
            derr = A[.,i] - Ax
            step = -(hg) * dx' :/ ((derr' * derr) + eta * (dx * dx'))
            conststep = min((1, max((0, step)))) 
            x = x + conststep * dx  
            err = (A, b) * (x' \ -1)
            vals[1,t] = zeta^2 * (x * x') + (err' * err) / row
            if (t>1) {
                dd = vals[1,t-1] - vals[1,t]
            }
        }
    }
    return(x)
}
end

*fw.step 
mata:
real vector fw(matrix A, matrix b, matrix x, eta) {    
    Ax = A * x'	
    hg = (Ax - b)' * A + eta * x
    i = select((1..cols(hg)), colmin(hg :== min(hg)))[1,1]
    dx = -x
    dx[1,i] = 1 - x[1,i]
    v = abs(min(dx))+abs(max(dx))
    if (v==0) {
        x = x
    }
    else {
        derr = A[.,i] - Ax
        step = -(hg) * dx' :/ ((derr' * derr) + eta * (dx * dx'))
        conststep = min((1, max((0, step)))) 
        x = x + conststep * dx  
    }
    return(x)
}
end

*spar function
mata:
    real matrix sspar(matrix V) {
        W = J(1,length(V),.)
        for (i=1; i<=length(V); ++i) {
            W[1,i] = V[1,i] <= max(V)/4 ? 0 : V[1,i]
        }
        W = W :/ sum(W)
        return(W)
    }
end
	
*sum normalize
mata:
    real vector sum_norm(matrix O) {
        sum_o = sum(O)
        if (sum_o!=0) {
            O = O / sum_o
        }
        else {
            O = J(1, cols(O), 1/cols(O)) 
        }
        return(O)
    }
end

*simple merge bt two vectors
mata:
    real matrix smerge(matrix A, matrix B) {
        v = J(rows(A), 1, .)
        A = (A, v)
        for (i=1; i<=rows(A); i++) {
            for (j=1; j<=rows(B); j++) {
                if (A[i,1]==B[j,1]) A[i,2]=B[j,2]
            }
	    }
        A = A[.,2]
        return(A)
}
end	

*projected covariates
mata:
    void projected(Y, Yprojected, Beta) {
        K = cols(Y)
        cdat = Y[selectindex(Y[,6]:==0),(1,2,4,8..K)]
        cdat = select(cdat, rowmissing(cdat):==0)
        X = cdat[.,4::cols(cdat)]
        NX = cols(X)
        yearFEs = uniqrows(cdat[.,3])
        for (fe=1;fe<=rows(yearFEs);fe++) {
            fevar = cdat[.,3]:==yearFEs[fe]
            X = (X,fevar)
        }
        unitFEs = uniqrows(cdat[.,2])
        for (fe=1;fe<=rows(unitFEs);fe++) {
            fevar = cdat[.,2]:==unitFEs[fe]
            X = (X,fevar)
        }            
        y = cdat[.,1]
        XX = quadcross(X,1 , X,1)
        Xy = quadcross(X,1 , y,0)
        beta = invsym(XX)*Xy
        Beta = beta[1::NX]
        X = Y[.,8::K]
        Yprojected = Y[.,1]-X*Beta
}
end

