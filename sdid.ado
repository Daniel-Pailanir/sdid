*! sdid: Synthetic Difference in Difference
*! Version 0.1.0 January 25, 2022
*! Author: Pailañir Daniel, Clarke Damian
*! dpailanir@fen.uchile.cl, dclarke@fen.uchile.cl

cap program drop sdid
program sdid
version 13.0
	
#delimit ;
    syntax varlist(min=4 numeric),
    [
    ]
    ;
#delimit cr  

*--------------------------------------------------------*
*- Create some temporal variables and locals            -*
*--------------------------------------------------------*
tokenize `varlist'
tempvar id diff tr
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

di as text "                          "
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

*di " ZetaOmega `ZetaOmega'"
*di " ZetaLambda `ZetaLambda'"

*-------------------------------------------------------*
*- Preparing data                                      -*
*-------------------------------------------------------*
*matrix of control units
preserve
qui keep if `tr'==0
keep `1' `id' `3'
qui levelsof `3', local(times) //local of all times
qui levelsof `3' if `3'<=`T0', local(timespre) //local of pre times
qui reshape wide `1', i(`id') j(`3')
mkmat _all, matrix(Y0) //matrix of control units Y0
restore

*matrix of treated units
preserve
qui keep if `tr'==1
keep `1' `id' `3'
qui reshape wide `1', i(`id') j(`3')
mkmat _all, matrix(Y1)
restore

preserve //for keep original data
matrix Y = (Y0 \ Y1) //matrix of control and treated units
clear
qui svmat Y
drop Y1
gen id = _n

local i=2
foreach n of local times {
    ren Y`i' t`n'
    local ++i
}

mkmat _all, matrix(Yall) //ver como hacer esto mejor, de guarda la bbdd 
mata : Yall = st_matrix("Yall")

egen promt = rowmean(t`Ttrmin'-t`T')
drop t`Ttrmin'-t`T' //drop post periods

local r=`N'+1
qui set obs `r'
  
forvalues t=`Tmin'/`T0' {
    qui sum     t`t' if id>`N0'
    qui replace t`t' = r(mean) in `r'
}

qui drop if id>`N0' & id!=. //drop treated units

mkmat _all, matrix(Y)
restore //end for keep original data

*-------------------------------------------------------*
*Matrices for optimization
*-------------------------------------------------------*
*Matrix A : Lambda
preserve
clear
qui svmat Y, names(col)

foreach t of local timespre {
    qui sum t`t' if id<=`N0'
    qui replace t`t' = t`t' - r(mean) 
}

drop id promt
qui keep in 1/`N0'
mkmat _all, matrix(A_l)
local col_l = colsof(A_l)
local row_l = rowsof(A_l)
mata : A_l = st_matrix("A_l")
restore

*Matrix b : Lambda
preserve
clear
qui svmat Y, names(col)
qui sum promt if id<=`N0'
qui replace promt = promt - r(mean)
keep promt
qui keep in 1/`N0'
mkmat _all, matrix(b_l)
mata : b_l = st_matrix("b_l")
restore

*Matrix A : Omega
preserve
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

mkmat _all, matrix(A_om)
keep t1-t`N0'
mkmat _all, matrix(A_o)
local col_o = colsof(A_o)
local row_o = rowsof(A_o)
mata : A_o = st_matrix("A_o")
restore

*Matrix b : Omega
preserve
clear
qui svmat A_om, names(col)
drop t1-t`N0'
mkmat _all, matrix(b_o)
mata : b_o = st_matrix("b_o")
restore

*eta value and omega vector
local eta_o = `row_o' * `ZetaOmega'^2
mata : lambda_o = J(1, `col_o', 1 / `col_o')

*eta value and lambda vector
local eta_l = `row_l' * `ZetaLambda'^2
mata : lambda_l = J(1, `col_l', 1 / `col_l')

*-------------------------------------------------------*
*LAMBDA
*-------------------------------------------------------*
local t=0
local maxIter=100
local mindecrease=(1e-5 * `sig')^2
mata : vals_l = J(1, `maxIter', .)
local dd=1

while (`t'<`maxIter' & (`t'<2 | `dd'>`mindecrease')) {
	local ++t    
    mata : Ax = A_l * lambda_l'	
    mata : hg = (Ax - b_l)' * A_l + `eta_l' * lambda_l
    mata : st_local("i", strofreal(select((1..cols(hg)), colmin(hg :== min(hg)))))
    mata :  dx = -lambda_l
    mata :  dx[1,`i'] = 1 - lambda_l[1,`i']
	mata : v=abs(min(dx))+abs(max(dx))
	mata : st_numscalar("v", v)

    if v==0 {
        mata : lambda_l = lambda_l
        mata : err = (A_l,b_l) * (lambda_l' \ -1)
        mata : vals_l[1, `t'] = `ZetaLambda'^2 * (lambda_l * lambda_l') + (err' * err) / `row_l'
    }
	else {
        mata : derr = A_l[1..`row_l',`i'] - Ax
        mata : step = -(hg) * dx' :/ ((derr' * derr) + `eta_l' * (dx * dx'))
        mata : st_local("step", strofreal(step))
        local conststep = min(1, max(0, `step'))
        mata : lambda_l = lambda_l + `conststep' * dx  
        mata : err = (A_l, b_l) * (lambda_l' \ -1)
    	mata : vals_l[1, `t'] = `ZetaLambda'^2 * (lambda_l * lambda_l') + (err' * err) / `row_l'
		
		if `t'>1 {
            mata : dd = vals_l[1, `t'-1] - vals_l[1, `t']
            mata : st_local("dd", strofreal(dd))
		}
    }
}

mata : st_local("maxlambda_l", strofreal(max(lambda_l)))
local cut = `maxlambda_l' / 4
mata : lambda_l = mm_cond(lambda_l :<=`cut', 0, lambda_l) //moremata install, creo
mata : lambda_l = lambda_l :/ sum(lambda_l)

local t=0
local maxIter=10000
local mindecrease=(1e-5 * `sig')^2
mata : vals_l = J(1, `maxIter', .)
local dd=1

while (`t'<`maxIter' & (`t'<2 | `dd'>`mindecrease')) {
	local ++t    
    mata : Ax = A_l * lambda_l'	
    mata : hg = (Ax - b_l)' * A_l + `eta_l' * lambda_l
    mata : st_local("i", strofreal(select((1..cols(hg)), colmin(hg :== min(hg)))))
    mata :  dx = -lambda_l
    mata :  dx[1,`i'] = 1 - lambda_l[1,`i']
	mata : v=abs(min(dx))+abs(max(dx))
	mata : st_numscalar("v", v)

    if v==0 {
        mata : lambda_l = lambda_l
        mata : err = (A_l,b_l) * (lambda_l' \ -1)
        mata : vals_l[1, `t'] = `ZetaLambda'^2 * (lambda_l * lambda_l') + (err' * err) / `row_l'
    }
	else {
        mata : derr = A_l[1..`row_l',`i'] - Ax
        mata : step = -(hg) * dx' :/ ((derr' * derr) + `eta_l' * (dx * dx'))
        mata : st_local("step", strofreal(step))
        local conststep = min(1, max(0, `step'))
        mata : lambda_l = lambda_l + `conststep' * dx  
        mata : err = (A_l, b_l) * (lambda_l' \ -1)
    	mata : vals_l[1, `t'] = `ZetaLambda'^2 * (lambda_l * lambda_l') + (err' * err) / `row_l'
		
		if `t'>1 {
            mata : dd = vals_l[1, `t'-1] - vals_l[1, `t']
            mata : st_local("dd", strofreal(dd))
		}
    }
}



*-------------------------------------------------------*
*OMEGA
*-------------------------------------------------------*
local t=0
local maxIter=100
local mindecrease=(1e-5 * `sig')^2
mata : vals_o = J(1, `maxIter', .)
local dd=1

while (`t'<`maxIter' & (`t'<2 | `dd'>`mindecrease')) {
	local ++t    
    mata : Ax = A_o * lambda_o'	
    mata : hg = (Ax - b_o)' * A_o + `eta_o' * lambda_o
    mata : st_local("i", strofreal(select((1..cols(hg)), colmin(hg :== min(hg)))))
    mata :  dx = -lambda_o
    mata :  dx[1,`i'] = 1 - lambda_o[1,`i']
	mata : v=abs(min(dx))+abs(max(dx))
	mata : st_numscalar("v", v)

    if v==0 {
        mata : lambda_o = lambda_o
        mata : err = (A_o,b_o) * (lambda_o' \ -1)
        mata : vals_o[1, `t'] = `ZetaOmega'^2 * (lambda_o * lambda_o') + (err' * err) / `row_o'
    }
	else {
        mata : derr = A_o[1..`row_o',`i'] - Ax
        mata : step = -(hg) * dx' :/ ((derr' * derr) + `eta_o' * (dx * dx'))
        mata : st_local("step", strofreal(step))
        local conststep = min(1, max(0, `step'))
        mata :  lambda_o = lambda_o + `conststep' * dx  
        mata :  err = (A_o, b_o) * (lambda_o' \ -1)
    	mata :  vals_o[1, `t'] = `ZetaOmega'^2 * (lambda_o * lambda_o') + (err' * err) / `row_o'

		if `t'>1 {
            mata : dd = vals_o[1, `t'-1] - vals_o[1, `t']
            mata : st_local("dd", strofreal(dd))
		}
		
    }
}

mata : st_local("maxlambda_o", strofreal(max(lambda_o)))
local cut = `maxlambda_o' / 4
mata : lambda_o = mm_cond(lambda_o :<=`cut', 0, lambda_o) //moremata install
mata : lambda_o = lambda_o :/ sum(lambda_o)

local t=0
local maxIter=10000
local mindecrease=(1e-5 * `sig')^2
mata : vals_o = J(1, `maxIter', .)
local dd=1

while (`t'<`maxIter' & (`t'<2 | `dd'>`mindecrease')) {
	local ++t    
    mata : Ax = A_o * lambda_o'	
    mata : hg = (Ax - b_o)' * A_o + `eta_o' * lambda_o
    mata : st_local("i", strofreal(select((1..cols(hg)), colmin(hg :== min(hg)))))
    mata :  dx = -lambda_o
    mata :  dx[1,`i'] = 1 - lambda_o[1,`i']
	mata : v=abs(min(dx))+abs(max(dx))
	mata : st_numscalar("v", v)

    if v==0 {
        mata : lambda_o = lambda_o
        mata : err = (A_o,b_o) * (lambda_o' \ -1)
        mata : vals_o[1, `t'] = `ZetaOmega'^2 * (lambda_o * lambda_o') + (err' * err) / `row_o'
    }
	else {
        mata : derr = A_o[1..`row_o',`i'] - Ax
        mata : step = -(hg) * dx' :/ ((derr' * derr) + `eta_o' * (dx * dx'))
        mata : st_local("step", strofreal(step))
        local conststep = min(1, max(0, `step'))
        mata :  lambda_o = lambda_o + `conststep' * dx  
        mata :  err = (A_o, b_o) * (lambda_o' \ -1)
    	mata :  vals_o[1, `t'] = `ZetaOmega'^2 * (lambda_o * lambda_o') + (err' * err) / `row_o'

		if `t'>1 {
            mata : dd = vals_o[1, `t'-1] - vals_o[1, `t']
            mata : st_local("dd", strofreal(dd))
		}
		
    }
}




*-------------------------------------------------------*
*TAU
*-------------------------------------------------------*
mata : tau = (-lambda_o, J(1, `Ntr', 1/`Ntr')) * Yall[1..`N',1..`Tobs'] * (-lambda_l, J(1, `Tpost', 1/`Tpost'))'

mata : st_local("tau", strofreal(tau))

di "{bf : objetivo post sparsify} : -15.60383"
di "{bf : tau post sparsify} : `tau'"

end

