*! sdid: Synthetic Difference in Difference
*! Version 0.1.0 January 25, 2022
*! Author: Pailañir Daniel, Clarke Damian
*! dpailanir@fen.uchile.cl, dclarke@fen.uchile.cl

cap program drop sdid
program sdid, eclass
version 13.0
	
#delimit ;
    syntax varlist(min=4 numeric), vce(string)
    [
    seed(integer 0)
    reps(integer 0)
    ]
    ;
#delimit cr  

*--------------------------------------------------------*
*- Create some temporal variables and locals            -*
*--------------------------------------------------------*
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
di as text "{c TLC}{hline 26}{c TT}{hline 6}{c TRC}"
di as text "{c |} Number of units          {c |} `N'   {c |}"
di as text "{c |} Number of times          {c |} `Tobs'   {c |}"
di as text "{c |} Smallest unit of time    {c |} `Tmin' {c |}"
di as text "{c |} Larger unit of time      {c |} `T' {c |}"
di as text "{c |} Number of treated units  {c |} `Ntr'    {c |}"
di as text "{c |} Number of control units  {c |} `N0'   {c |}"
di as text "{c |} Number of post-periods   {c |} `Tpost'   {c |}"
di as text "{c |} Number of pre-periods    {c |} `Tpre'   {c |}"
di as text "{c |} Maximun time of control  {c |} `T0' {c |}"
di as text "{c |} First time of treatment  {c |} `Ttrmin' {c |}"
di as text "{c BLC}{hline 26}{c BT}{hline 6}{c BRC}"

*-------------------------------------------------------*
*- Calculate \zeta                                     -*
*-------------------------------------------------------*
bys `id' : egen `tr' = mean(`4')
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

*restore original data
use `data', clear

*------------------------------------------------------------------------------*
*VCE : bootstrap
*------------------------------------------------------------------------------*
if "`vce'"=="bootstrap" {
    set seed `seed'
    local b = 1
    local B = `reps'
    mata: tau_b = J(1, `B', .)

    if (`Ntr'==1)==1 {
        di as err "It is not possible to do Bootstrap se because there is only one treated unit"
        exit 198
    }

    *Sum normalize for initialization
    mata: lambda_o = sum_norm(lambda_o)
		
while `b'<=`B' {
    preserve
    bsample , cluster(`id') idcluster(`id2')
    bys `id2' : egen `tr'`b' = mean(`4')
    qui replace `tr'`b' = 1 if `tr'`b'!=0

    qui sum `tr'`b'
    if (r(mean)==0 | r(mean)==1) {
        *di "Boot `b' : Nothing to do"
    }
    else {
        *di "Boot `b' : Running"
        *-------------------------------------------------------*
        *- Preparing data                                      -*
        *-------------------------------------------------------*
        *original bootstrap data
        tempfile data`b'
        qui save "`data`b''"
        qui levelsof `3', local(times)
        qui levelsof `3' if `3'<=`T0', local(timespre)

        *matrix of control units
        qui keep if `tr'`b'==0
        keep `1' `id2' `3'
        qui reshape wide `1', i(`id2') j(`3')
        mkmat _all, matrix(Y0)

        *matrix of treated units
        use `data`b'', clear
        qui keep if `tr'`b'==1
        keep `1' `id2' `3'
        qui reshape wide `1', i(`id2') j(`3')
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

        *data for estimator
        mkmat _all, matrix(Yall`b')
        mata: Yall`b' = st_matrix("Yall`b'")

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
        *Matrices for optimization
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
        mata : A_o = st_matrix("A_o")
        mata : b_o = st_matrix("b_o")

        *eta value
        local eta_o = `row_o' * `ZetaOmega'^2
        local eta_l = `row_l' * `ZetaLambda'^2
        *----------------------------------------------------------------------*
        *LAMBDA
        *----------------------------------------------------------------------*
        local mindec = (1e-5 * `sig')^2
        mata: w_l = lambda(A_l, b_l, lambda_l, `eta_l', `ZetaLambda', 100, `mindec')
        mata: w_l = sspar(w_l)
        mata: w_l = lambda(A_l, b_l, w_l, `eta_l', `ZetaLambda', 10000, `mindec')
        *----------------------------------------------------------------------*
        *OMEGA
        *-----------------------------------------------------------------------*
        mata: w_o = lambda(A_o, b_o, lambda_o, `eta_o', `ZetaOmega', 100, `mindec')
        mata: w_o = sspar(w_o)
        mata: w_o = lambda(A_o, b_o, w_o, `eta_o', `ZetaOmega', 10000, `mindec')
        *-----------------------------------------------------------------------*
        *TAU
        *-----------------------------------------------------------------------*
        mata: tau_b[1,`b'] = (-w_o, J(1, `Ntr', 1/`Ntr'))*Yall`b'[1..`N',1..`Tobs']*(-w_l, J(1, `Tpost', 1/`Tpost'))'
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
    
    *Sum normalize for initialization
    mata: lambda_o = sum_norm(lambda_o)
		
    while `b'<=`B' {
        preserve
        clear
        tempvar i ind
        qui set obs `N0'
		gen `i' = _n
        gen `ind' = runiform()
        sort `ind'
        local nN0 = `N0' - `Ntr'
        qui putmata ind1 = `i', replace
        qui drop in `N0'
        qui levelsof `i', local(in)	
        mata: ind2 = ind1[1..`nN0',1]

        use `data', clear
        qui xtset `id' `3'
        qui drop if `id'>`N0'
        local sum1 = 0
        foreach s of local in {
            local sum1 = `sum1' + `s' 
        }

        local sum2 = 0
        qui levelsof `id', local(j)
        foreach s of local j {
            local sum2 = `sum2' + `s' 
        }

		local drp = `sum2'-`sum1'
        qui drop if `id'==`drp'
		
        *eta value
        mata: st_local("row_o", strofreal(rows(A_o[1..`Tpre',(ind2)])))
        mata: st_local("row_l", strofreal(rows(A_l[(ind2),1..`Tpre'])))
        local eta_o = `row_o' * `ZetaOmega'^2
        local eta_l = `row_l' * `ZetaLambda'^2
        *----------------------------------------------------------------------*
        *LAMBDA
        *----------------------------------------------------------------------*
        local mindec = (1e-5 * `sig')^2
        mata: w_l = lambda(A_l[(ind2),1..`Tpre'], b_l[(ind2),1], lambda_l, `eta_l', `ZetaLambda', 100, `mindec')
        mata: w_l = sspar(w_l)
        mata: w_l = lambda(A_l[(ind2),1..`Tpre'], b_l[(ind2),1], w_l, `eta_l', `ZetaLambda', 10000, `mindec')
        *----------------------------------------------------------------------*
        *OMEGA
        *----------------------------------------------------------------------*
        mata: w_o = lambda(A_o[1..`Tpre',(ind2)], b_o, lambda_o[1,(ind2)], `eta_o', `ZetaOmega', 100, `mindec')
        mata: w_o = sspar(w_o)
        mata: w_o = lambda(A_o[1..`Tpre',(ind2)], b_o, w_o, `eta_o', `ZetaOmega', 10000, `mindec')
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
if "`vce'"=="jackknife" {
    if (`Ntr'==1)==1 {
        di as err "It is not possible to do Jackknife se because there is only one treated unit"
        exit 198
    }
    
}


*Display results and save results
ereturn local se `se' 
ereturn local tau `tau'

di as text " "
di as text "{c TLC}{hline 16}{c TT}{hline 11}{c TRC}"
di as text "{c |} {bf: tau}           {c |} " as result %9.5f `tau'  as text " {c |}"
di as text "{c |} {bf: se bootstrap}  {c |} " as result %9.5f `se' as text " {c |}"
di as text "{c BLC}{hline 16}{c BT}{hline 11}{c BRC}"

*Restore original data
use `data', clear

end

mata:
function lambda(matrix A, matrix b, matrix x, eta, zeta, maxIter, mindecrease)
{

row = rows(A)
col = cols(A)
vals = J(1, maxIter, .)
t=0
dd=1

while (t<maxIter & (t<2 | dd>mindecrease)) {
    t=t+1    
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
        derr = A[1..row,i] - Ax
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

mata:
function sspar(matrix V)
{
    W = mm_cond(V :<= max(V)/4, 0, V)
    W = W :/ sum(W)
    return(W)
}
end
	
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
			
