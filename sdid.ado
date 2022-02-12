*! sdid: Synthetic Difference in Difference
*! Version 0.1.0 January 25, 2022
*! Author: Pailañir Daniel, Clarke Damian
*! dpailanir@fen.uchile.cl, dclarke@fen.uchile.cl

cap program drop sdid
program sdid, eclass
version 13.0
	
#delimit ;
    syntax varlist(min=4 numeric), vce(string) adoption(string)
    [
    seed(integer 0)
    reps(integer 0)
    ]
    ;
#delimit cr  

*------------------------------------------------------------------------------*
*ONE TIME ADOPTION
*------------------------------------------------------------------------------*
if "`adoption'"=="normal" {
tokenize `varlist'
tempvar id id2 diff tr
egen `id' = group(`2')

qui xtset `id' `3'
local N    = r(imax) //number of units
local Tmin = r(tmin) //t min
local T    = r(tmax) //t max

qui tab `id' if `4'==1
local Ntr = r(r)         //number of treated units
local N0  = `N' - `Ntr'  //number of control units

qui tab `3' if `4'==1
local Tpost  = r(r)             //number of post times
local T0     = `T'  - `Tpost'   //max time of control
local Tobs   = `T'  - `Tmin' +1 //number of times
local Tpre   = `T0' - `Tmin' +1 //number of pre times
local Ttrmin = `T0' + 1         //first year of treatment

di as text " "
di as text "{c TLC}{hline 26}{c TT}{hline 6}"
di as text "{c |} Number of units          {c |} `N'   "
di as text "{c |} Number of times          {c |} `Tobs'   "
di as text "{c |} Smallest unit of time    {c |} `Tmin' "
di as text "{c |} Larger unit of time      {c |} `T' "
di as text "{c |} Number of treated units  {c |} `Ntr'   "
di as text "{c |} Number of control units  {c |} `N0'   "
di as text "{c |} Number of post-periods   {c |} `Tpost'   "
di as text "{c |} Number of pre-periods    {c |} `Tpre'   "
di as text "{c |} Maximun time of control  {c |} `T0' "
di as text "{c |} First time of treatment  {c |} `Ttrmin' "
di as text "{c BLC}{hline 26}{c BT}{hline 6}"

*-------------------------------------------------------*
*- Calculate \zeta                                     -*
*-------------------------------------------------------*
bys `id': egen `tr' = mean(`4')
qui replace `tr' = 1 if `tr'!=0
local EtaOmega  = (`Ntr' * `Tpost')^(1/4)
local EtaLambda = 1e-6
qui gen `diff' = `1' - L.`1'
qui sum `diff' if `3'<=`T0' & `tr'==0
local sig = r(sd)
local ZetaOmega  = `EtaOmega'  * `sig' 
local ZetaLambda = `EtaLambda' * `sig'
*-------------------------------------------------------*
*- Preparing data                                      -*
*-------------------------------------------------------*
*original data
tempfile data
qui save "`data'"
qui levelsof `3', local(times)                 //local of all times
qui levelsof `3' if `3'<=`T0', local(timespre) //local of pre times

*matrix of control units
qui keep if `tr'==0
keep `1' `id' `3'
qui reshape wide `1', i(`id') j(`3')
mkmat _all, matrix(Y0)

*matrix of treated units
use `data', clear
qui keep if `tr'==1
keep `1' `id' `3'
qui reshape wide `1', i(`id') j(`3')
mkmat _all, matrix(Y1)

*matrix of control and treated units
matrix Y = (Y0 \ Y1) 
clear
qui svmat Y
gen id = _n
qui putmata idd1 = Y1 idd2 = id, replace
mata: indd = (idd1, idd2)
drop Y1

local i=2
foreach n of local times {
    ren Y`i' t`n'
    local ++i
}

*All data for estimator
mkmat _all, matrix(Yall)
mata: Yall = st_matrix("Yall")

egen promt = rowmean(t`Ttrmin'-t`T')
drop t`Ttrmin'-t`T'
local r=`N'+1
qui set obs `r'
  
forvalues t=`Tmin'/`T0' {
    qui sum     t`t' if id>`N0'
    qui replace t`t' = r(mean) in `r'
}

qui drop if id>`N0' & id!=.
mkmat _all, matrix(Y)
*-------------------------------------------------------*
*- Matrices for optimization                           -*
*-------------------------------------------------------*
*Matrix A & b : Lambda 
clear
qui svmat Y, names(col)

local vr `timespre' promt
foreach t of local vr {
    if "`t'"=="promt" local n ""
    if "`t'"!="promt" local n "t"
    qui sum `n'`t' if id<=`N0'
    qui replace `n'`t' = `n'`t' - r(mean) 
}

qui keep in 1/`N0'
mkmat promt, matrix(b_l)
mkmat t`Tmin'-t`T0', matrix(A_l)
local col_l = colsof(A_l)
local row_l = rowsof(A_l)
mata: A_l = st_matrix("A_l")
mata: b_l = st_matrix("b_l")

*Matrix A & b : Omega
clear
qui svmat Y, names(col)
drop promt id
gen id = _n
qui reshape long t, i(id) j(a)
qui reshape wide t, i(a) j(id)
drop a

local max=`N0'+1
forvalues t=1/`max'  {
    qui sum t`t'
    qui replace t`t' = t`t' - r(mean)
}

mkmat t`max', matrix(b_o)
mkmat t1-t`N0', matrix(A_o)
local col_o = colsof(A_o)
local row_o = rowsof(A_o)
mata: A_o = st_matrix("A_o")
mata: b_o = st_matrix("b_o")

*eta value and omega vector
local eta_o = `row_o' * `ZetaOmega'^2
local eta_l = `row_l' * `ZetaLambda'^2

*eta value and lambda vector
mata: lambda_o = J(1, `col_o', 1 / `col_o')
mata: lambda_l = J(1, `col_l', 1 / `col_l')
*------------------------------------------------------------------------------*
*LAMBDA
*------------------------------------------------------------------------------*
local mindecrease=(1e-5 * `sig')^2
mata: lambda_l = lambda(A_l, b_l, lambda_l, `eta_l', `ZetaLambda', 100, `mindecrease')
mata: lambda_l = sspar(lambda_l)
mata: lambda_l = lambda(A_l, b_l, lambda_l, `eta_l', `ZetaLambda', 10000, `mindecrease')
*------------------------------------------------------------------------------*
*OMEGA
*------------------------------------------------------------------------------*
mata: lambda_o = lambda(A_o, b_o, lambda_o, `eta_o', `ZetaOmega', 100, `mindecrease')
mata: lambda_o = sspar(lambda_o)
mata: lambda_o = lambda(A_o, b_o, lambda_o, `eta_o', `ZetaOmega', 10000, `mindecrease')

*save weights in e(r)
mata: st_matrix("lambda", lambda_l')
mata: st_matrix("omega", lambda_o')
ereturn matrix lambda lambda
ereturn matrix omega  omega
*------------------------------------------------------------------------------*
*TAU
*------------------------------------------------------------------------------*
mata: tau = (-lambda_o, J(1, `Ntr', 1/`Ntr'))*Yall[1..`N',1..`Tobs']*(-lambda_l, J(1, `Tpost', 1/`Tpost'))'
mata: st_local("tau", strofreal(tau))
*------------------------------------------------------------------------------*
*VCE : bootstrap
*------------------------------------------------------------------------------*
if "`vce'"=="bootstrap" {
    set seed `seed'
    local b = 1
    local B = `reps'
    mata: tau_b = J(`B', 1, .)

    if (`Ntr'==1)==1 {
        di as err "It is not possible to do Bootstrap se because there is only one treated unit"
        exit 198
    }
		
    while `b'<=`B' {
        preserve
        clear
        tempvar i1 i2
        qui set obs `N'
        gen `i1' = _n
        bsample , cluster(`i1') idcluster(`i2')
        qui tab `i1' if `i1'<=`N0'
        local r1 = r(r)
        qui tab `i1' if `i1'>`N0'
        local r2 = r(r)

        if (`r1'==0 | `r2'==0) {
            *di "all control or treated unit"
        }	
        else {
            *di "at least one treated and control unit"
            sort `i1'
            qui putmata ind1 = `i1', replace
            mata: ind2 = select(ind1, ind1:<=`N0')
            mata: Ntr = length(ind1) - length(ind2)
            mata: l_o = sum_norm(lambda_o[,ind2]) 

            *eta value
            mata: st_local("row_l", strofreal(rows(A_l[(ind2),1..`Tpre'])))
            local eta_o = `row_o' * `ZetaOmega'^2
            local eta_l = `row_l' * `ZetaLambda'^2
            *------------------------------------------------------------------*
            *LAMBDA
            *------------------------------------------------------------------*
            local mindec = (1e-5 * `sig')^2
            mata: w_l = lambda(A_l[(ind2),.], b_l[(ind2),1], lambda_l, `eta_l', `ZetaLambda', 100, `mindec')
            mata: w_l = sspar(w_l)
            mata: w_l = lambda(A_l[(ind2),.], b_l[(ind2),1], w_l, `eta_l', `ZetaLambda', 10000, `mindec')
            *------------------------------------------------------------------*
            *OMEGA
            *------------------------------------------------------------------*
            mata: w_o = lambda(A_o[.,(ind2)], b_o, l_o, `eta_o', `ZetaOmega', 100, `mindec')
            mata: w_o = sspar(w_o)
            mata: w_o = lambda(A_o[.,(ind2)], b_o, w_o, `eta_o', `ZetaOmega', 10000, `mindec')
            *------------------------------------------------------------------*
            *TAU
            *------------------------------------------------------------------*
            mata: tau_b[`b',] = (-w_o, J(1, Ntr, 1/Ntr))*Yall[(ind1),1..`Tobs']*(-w_l, J(1, `Tpost', 1/`Tpost'))'
           local ++b
        }
        restore
    }
    mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(tau_b)))
    mata: st_local("se", strofreal(se))
}

*------------------------------------------------------------------------------*
*VCE : placebo
*------------------------------------------------------------------------------*
else if "`vce'"=="placebo" {
    set seed `seed'
    local b = 1
    local B = `reps'
    mata : tau_p = J(1, `B', .)

    if (`N0'<=`Ntr')==1 {
        di as err "It is not possible to do Placebo se because there are more treated units than control units"
        exit 198
    }
    while `b'<=`B' {
        preserve
        clear
        tempvar i iu
        qui set obs `N0'
		gen `i' = _n
        gen `iu' = runiform()
        sort `iu'
		qui putmata ind1 = `i', replace
        local nN0 = `N0' - `Ntr'
        local mN0 = `nN0' + 1
        qui drop in `mN0'/`N0'
        mata: ind2 = ind1[1..`nN0',1]	
        mata: l_o = sum_norm(lambda_o[,ind2]) 

        *eta value
        mata: st_local("row_l", strofreal(rows(A_l[(ind2),1..`Tpre'])))
        local eta_o = `row_o' * `ZetaOmega'^2
        local eta_l = `row_l' * `ZetaLambda'^2
        *----------------------------------------------------------------------*
        *LAMBDA
        *----------------------------------------------------------------------*
        local mindec = (1e-5 * `sig')^2
        mata: w_l = lambda(A_l[(ind2),.], b_l[(ind2),1], lambda_l, `eta_l', `ZetaLambda', 100, `mindec')
        mata: w_l = sspar(w_l)
        mata: w_l = lambda(A_l[(ind2),.], b_l[(ind2),1], w_l, `eta_l', `ZetaLambda', 10000, `mindec')
        *----------------------------------------------------------------------*
        *OMEGA
        *----------------------------------------------------------------------*
        mata: w_o = lambda(A_o[.,(ind2)], b_o, l_o, `eta_o', `ZetaOmega', 100, `mindec')
        mata: w_o = sspar(w_o)
        mata: w_o = lambda(A_o[.,(ind2)], b_o, w_o, `eta_o', `ZetaOmega', 10000, `mindec')
        *----------------------------------------------------------------------*
        *TAU
        *----------------------------------------------------------------------*
        mata: tau_p[1,`b'] = (-w_o, J(1, `Ntr', 1/`Ntr'))*Yall[(ind1),1..`Tobs']*(-w_l, J(1, `Tpost', 1/`Tpost'))'
        local ++b
        restore
    }
    mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(tau_p)))
    mata: st_local("se", strofreal(se))
}

*------------------------------------------------------------------------------*
*VCE : jackknife
*------------------------------------------------------------------------------*
else if "`vce'"=="jackknife" {
    if (`Ntr'==1)==1 {
        di as err "It is not possible to do Jackknife se because there is only one treated unit"
        exit 198
    }
    mata: se = jackknife(Yall, lambda_l, lambda_o, `N0', `Tpost', `N', `Tobs')
    mata: st_local("se", strofreal(se))
}

*------------------------------------------------------------------------------*
*Display results and save results
*------------------------------------------------------------------------------*
ereturn local se `se' 
ereturn local tau `tau'

di as text " "
di as text "{c TLC}{hline 7}{c TT}{hline 11}{c TRC}"
di as text "{c |} {bf: tau}  {c |} " as result %9.5f `tau'  as text " {c |}"
di as text "{c |} {bf: se}   {c |} " as result %9.5f `se' as text " {c |}"
di as text "{c BLC}{hline 7}{c BT}{hline 11}{c BRC}"
}

*------------------------------------------------------------------------------*
*STAGGERED ADOPTION
*------------------------------------------------------------------------------*
else if "`adoption'"=="staggered" {
    tokenize `varlist'
    tempvar i m adoption trt
    egen `i' = group(`2')
    qui xtset `i' `3'
    local Tmin = r(tmin) //t min
    local T    = r(tmax) //t max
    *Define treated periods
    bys `2' `4': egen `m' = min(`3')
    by  `2': egen `adoption' = max(`m')
    qui by  `2': replace `adoption' = 0 if `adoption'==`Tmin'

    qui levelsof `adoption' if `adoption'>0, local(trt)
    qui tab `adoption' if `adoption'>0
    local lenght = r(r) 
    mata: results = J(`lenght', 6, .)

*filter data and SDiD
mata: i = 1
foreach t of local trt {
    preserve
    tempvar id id2 diff tr post_treat
    qui keep if `adoption'==`t' | `adoption'==0
    qui sum `adoption'
    scalar time=r(max)
    gen `post_treat' = 0
    qui replace `post_treat' = 1 if `3'>=`adoption' & `adoption'!=0
    qui sum `post_treat' if `post_treat'==1
    scalar obs=r(N)
    *--------------------------------------------------------*
    *- Create some temporal variables and locals            -*
    *--------------------------------------------------------*
    egen `id' = group(`2')
    qui xtset `id' `3'
    local N    = r(imax)      //number of units
    qui tab `id' if `4'==1
    local Ntr = r(r)         //number of treated units
    local N0  = `N' - `Ntr'  //number of control units
    qui tab `3' if `3'>=`t'
    local Tpost  = r(r)             //number of post times
    local T0     = `T'  - `Tpost'   //max time of control
    local Tobs   = `T'  - `Tmin' +1 //number of times
    local Tpre   = `T0' - `Tmin' +1 //number of pre times
    local Ttrmin = `T0' + 1         //first year of treatment
    *-------------------------------------------------------*
    *- Calculate \zeta                                     -*
    *-------------------------------------------------------*
    bys `id': egen `tr' = mean(`4')
    qui replace `tr' = 1 if `tr'!=0
    local EtaOmega  = (`Ntr' * `Tpost')^(1/4)
    local EtaLambda = 1e-6
    qui gen `diff' = `1' - L.`1'
    qui sum `diff' if `3'<=`T0' & `tr'==0
    local sig = r(sd)
    local ZetaOmega  = `EtaOmega'  * `sig' 
    local ZetaLambda = `EtaLambda' * `sig'
    *-------------------------------------------------------*
    *- Preparing data                                      -*
    *-------------------------------------------------------*
    *original data
    tempfile data
    qui save "`data'"
    qui levelsof `3', local(times)                 
    qui levelsof `3' if `3'<=`T0', local(timespre) 
    *matrix of control units
    qui keep if `tr'==0
    keep `1' `id' `3'
    qui reshape wide `1', i(`id') j(`3')
    mkmat _all, matrix(Y0)
    *matrix of treated units
    use `data', clear
    qui keep if `tr'==1
    keep `1' `id' `3'
    qui reshape wide `1', i(`id') j(`3')
    mkmat _all, matrix(Y1)
    *matrix of control and treated units
    matrix Y = (Y0 \ Y1) 
    clear
    qui svmat Y
    drop Y1
    gen id = _n
    local i=2
    foreach n of local times {
        ren Y`i' t`n'
        local ++i
    }

    *All data for estimator
    mkmat _all, matrix(Yall)
    mata: Yall = st_matrix("Yall")
    egen promt = rowmean(t`Ttrmin'-t`T')
    drop t`Ttrmin'-t`T'
    local r=`N'+1
    qui set obs `r'
  
    forvalues t=`Tmin'/`T0' {
        qui sum     t`t' if id>`N0'
        qui replace t`t' = r(mean) in `r'
    }

    qui drop if id>`N0' & id!=.
    mkmat _all, matrix(Y)
    *-------------------------------------------------------*
    *- Matrices for optimization                           -*
    *-------------------------------------------------------*
    *Matrix A & b : Lambda 
    clear
    qui svmat Y, names(col)

    local vr `timespre' promt
    foreach t of local vr {
        if "`t'"=="promt" local n ""
        if "`t'"!="promt" local n "t"
        qui sum `n'`t' if id<=`N0'
        qui replace `n'`t' = `n'`t' - r(mean) 
    }

    qui keep in 1/`N0'
    mkmat promt, matrix(b_l)
    mkmat t`Tmin'-t`T0', matrix(A_l)
    local col_l = colsof(A_l)
    local row_l = rowsof(A_l)
    mata: A_l = st_matrix("A_l")
    mata: b_l = st_matrix("b_l")
    *Matrix A & b : Omega
    clear
    qui svmat Y, names(col)
    drop promt id
    gen id = _n
    qui reshape long t, i(id) j(a)
    qui reshape wide t, i(a) j(id)
    drop a

    local max=`N0'+1
    forvalues t=1/`max'  {
        qui sum t`t'
        qui replace t`t' = t`t' - r(mean)
    }
    mkmat t`max', matrix(b_o)
    mkmat t1-t`N0', matrix(A_o)
    local col_o = colsof(A_o)
    local row_o = rowsof(A_o)
    mata: A_o = st_matrix("A_o")
    mata: b_o = st_matrix("b_o")

    *eta value and omega vector
    local eta_o = `row_o' * `ZetaOmega'^2
    local eta_l = `row_l' * `ZetaLambda'^2

    *eta value and lambda vector
    mata: lambda_o = J(1, `col_o', 1 / `col_o')
    mata: lambda_l = J(1, `col_l', 1 / `col_l')		
    *--------------------------------------------------------------------------*
    *LAMBDA
    *--------------------------------------------------------------------------*
    local mindecrease=(1e-5 * `sig')^2
    mata: lambda_l = lambda(A_l, b_l, lambda_l, `eta_l', `ZetaLambda', 100, `mindecrease')
    mata: lambda_l = sspar(lambda_l)
    mata: lambda_l = lambda(A_l, b_l, lambda_l, `eta_l', `ZetaLambda', 10000, `mindecrease')
    *------------------------------------------------------------------------------*
    *OMEGA
    *------------------------------------------------------------------------------*
    mata: lambda_o = lambda(A_o, b_o, lambda_o, `eta_o', `ZetaOmega', 100, `mindecrease')
    mata: lambda_o = sspar(lambda_o)
    mata: lambda_o = lambda(A_o, b_o, lambda_o, `eta_o', `ZetaOmega', 10000, `mindecrease')
    *save weights in e(r)
    mata: st_matrix("lambda", lambda_l')
    mata: st_matrix("omega", lambda_o')
    ereturn matrix lambda lambda
    ereturn matrix omega  omega
    *--------------------------------------------------------------------------*
    *TAU
    *--------------------------------------------------------------------------*
    mata: results[i,1] = st_numscalar("time")
    mata: results[i,2] = st_numscalar("obs")
    mata: results[i,3] = (-lambda_o, J(1, `Ntr', 1/`Ntr'))*Yall[1..`N',1..`Tobs']*(-lambda_l, J(1, `Tpost', 1/`Tpost'))'
    mata: i = i+1
    restore
}
    scalar drop time obs
    mata: results[.,2] = results[.,2]/sum(results[.,2])
    mata: ATT = sum(results[.,2] :* results[.,3])
    mata: st_local("ATT", strofreal(ATT))
	
*Display results
di as text " "
di as text "{c TLC}{hline 8}{c TT}{hline 11}{c TRC}"
di as text "{c |} {bf: ATT}   {c |} " as result %9.5f `ATT'  as text " {c |}"
di as text "{c |} {bf: se}    {c |} " as text "          {c |}"
di as text "{c BLC}{hline 8}{c BT}{hline 11}{c BRC}"
	
}


*Restore original data
use `data', clear

end

*minimization
mata:
function lambda(matrix A, matrix b, matrix x, eta, zeta, maxIter, mindecrease)
{

row = rows(A)
col = cols(A)
vals = J(1, maxIter, .)
t=0
dd=1

while (t<maxIter & (t<2 | dd>mindecrease)) {
    t++
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

*spar function
mata:
function sspar(matrix V)
{
    W = mm_cond(V :<= max(V)/4, 0, V)
    W = W :/ sum(W)
    return(W)
}
end
	
*sum normalize
mata:
function sum_norm(matrix O)
{
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
function smerge(matrix A, matrix B)
{
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

*for jackknife
mata:
function jackknife(matrix Y, matrix L, matrix O, c, tp, N, T)
{
    tau_j = J(N, 1, .)
	ind1 = (1::N)
    for (i=1; i<=N; i++) {
        id = select(ind1, ind1:!=i)
        id2 = select(id, id:<=c)
        t = length(id)-length(id2)
		l_o = sum_norm(O[,id2]) 
        tau_j[i,1] = (-l_o, J(1, t, 1/t))*Y[(id),1..T]*(-L, J(1, tp, 1/tp))'
    }
    se_j = sqrt(((N-1)/N) * (N - 1) * variance(vec(tau_j)))
    return(se_j)
}
end
		
			
