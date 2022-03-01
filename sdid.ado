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
    reps(integer 50)
    ]
    ;
#delimit cr  

/*
To do:
 (1) Error check (ensure blance, missings, ...)
 (2) Standard errors for staggered adoption
 (3) Graphing of results 
 (4) Incorporate controls
*/

*------------------------------------------------------------------------------*
* (0) Error checks in parsing
*------------------------------------------------------------------------------*


preserve
*------------------------------------------------------------------------------*
* (1) ONE TIME ADOPTION
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
    di as text "{c |} Maximum time of control  {c |} `T0' "
    di as text "{c |} First time of treatment  {c |} `Ttrmin' "
    di as text "{c BLC}{hline 26}{c BT}{hline 6}"
    
    *--------------------------------------------------------------------------*
    *- Calculate \zeta                                                        -*
    *--------------------------------------------------------------------------*
    bys `id': egen `tr' = mean(`4')
    qui replace `tr' = 1 if `tr'!=0
    local EtaOmega  = (`Ntr' * `Tpost')^(1/4)
    local EtaLambda = 1.0c6f7a0b5ed8dX-014
    qui gen `diff' = `1' - L.`1'
    qui sum `diff' if `3'<=`T0' & `tr'==0
    local sig = r(sd)
    local ZetaOmega  = `EtaOmega'  * `sig' 
    local ZetaLambda = `EtaLambda' * `sig'
    
    *--------------------------------------------------------------------------*
    *- Preparing data                                                         -*
    *--------------------------------------------------------------------------*
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
    mata: Yall = Yall[1..`N',1..`Tobs']
    
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
    *--------------------------------------------------------------------------*
    *- Matrices for optimization                                              -*
    *--------------------------------------------------------------------------*
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
    
    local mindec=(1.4f8b588e368f1X-011 * `sig')^2
    local col_l = colsof(A_l)
    local col_o = colsof(A_o)
    mat def lambda_l = J(1, `col_l', 1 / `col_l')
    mat def lambda_o = J(1, `col_o', 1 / `col_o')

    mata: A_l=st_matrix("A_l")
    mata: b_l=st_matrix("b_l")
    mata: A_o=st_matrix("A_o")
    mata: b_o=st_matrix("b_o")
    
    #delimit ;
    mata: tau = estTau(A_l, b_l, A_o, b_o, st_matrix("lambda_l"),st_matrix("lambda_o"),
                       Yall, `ZetaLambda', `ZetaOmega',`mindec',`Ntr',`N',
                       `Tobs',`Tpost');
    #delimit cr
    ereturn matrix lambda lambda
    ereturn matrix omega  omega
    matrix LAMBDA_L = e(lambda)
    matrix LAMBDA_O = e(omega)
    mata: st_local("tau", strofreal(tau))
    
    *--------------------------------------------------------------------------*
    *Standard error: bootstrap
    *--------------------------------------------------------------------------*
    if "`vce'"=="bootstrap" {
        set seed `seed'
        local b = 1
        local B = `reps'
        mata: tau_b = J(`B', 1, .)
    
        if (`Ntr'==1)==1 {
            di as err "It is not possible to do Bootstrap se because there is only one treated unit"
            exit 198
        }
        dis "Bootstrap replications (`reps'). This may take some time."
        dis "----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5"

        clear
        while `b'<=`B' {
            tempvar i1 i2
            qui set obs `N'
            gen `i1' = _n
            bsample, cluster(`i1') idcluster(`i2')
            qui tab `i1' if `i1'<=`N0'
            local r1 = r(r)
            qui tab `i1' if `i1'>`N0'
            local r2 = r(r)
    
            if (`r1'==0 | `r2'==0) {
                *di "all control or treated unit"
            }	
            else {
                display in smcl "." _continue
                if mod(`b',50)==0 dis "     `b'"
                
                sort `i1'
                qui putmata ind1 = `i1', replace
                mata: ind2 = select(ind1, ind1:<=`N0')
                mata: Ntrb = length(ind1) - length(ind2)
                mata: lambda_o=st_matrix("LAMBDA_O")
                mata: l_o = sum_norm(lambda_o'[,ind2]) 
                mata: A_lb = A_l[(ind2),.]
                mata: b_lb = b_l[(ind2),1] 
                mata: A_ob = A_o[.,(ind2)]
                mata: Yallb = Yall[(ind1),1..`Tobs']
                
                local mindec = (1e-5 * `sig')^2
                #delimit ;
                mata: tau_b[`b',] = estTau(A_lb, b_lb, A_ob, b_o, st_matrix("LAMBDA_L")',
                                           l_o, Yallb, `ZetaLambda', `ZetaOmega',
                                           `mindec',Ntrb,`N',`Tobs',`Tpost');
                #delimit cr                
                local ++b
            }
        
        }
        mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(tau_b)))
        mata: st_local("se", strofreal(se))
    }
    
    *--------------------------------------------------------------------------*
    *Standard error: placebo
    *--------------------------------------------------------------------------*
    else if "`vce'"=="placebo" {
        set seed `seed'
        local b = 1
        local B = `reps'
        mata : tau_p = J(1, `B', .)
    
        if (`N0'<=`Ntr')==1 {
            di as err "It is not possible to do Placebo se because there are more treated units than control units"
            exit 198
        }
        else {
            dis "Placebo inference (`reps' permutations). This may take some time."
            dis "----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5"
        }
        while `b'<=`B' {
            display in smcl "." _continue
            if mod(`b',50)==0 dis "     `b'"
            
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
            mata: lambda_o=st_matrix("LAMBDA_O")'
            mata: l_o = sum_norm(lambda_o[,ind2])
            
            mata: A_lb = A_l[(ind2),.]
            mata: b_lb = b_l[(ind2),1] 
            mata: A_ob = A_o[.,(ind2)]
            mata: Yallb = Yall[(ind1),1..`Tobs']

            local mindec = (1e-5 * `sig')^2
            #delimit ;
            mata: tau_p[1,`b'] = estTau(A_lb, b_lb, A_ob, b_o, st_matrix("LAMBDA_L")',
                                        l_o, Yallb, `ZetaLambda', `ZetaOmega',
                                        `mindec',`Ntr',`N',`Tobs',`Tpost');
            #delimit cr
            local ++b
        }
        mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(tau_p)))
        mata: st_local("se", strofreal(se))
    }
    
    *--------------------------------------------------------------------------*
    *Standard error: jackknife
    *--------------------------------------------------------------------------*
    else if "`vce'"=="jackknife" {
        if (`Ntr'==1)==1 {
            di as err "It is not possible to do Jackknife se because there is only one treated unit"
            exit 198
        }
        mata: se = jackknife(Yall, st_matrix("LAMBDA_L")', st_matrix("LAMBDA_O")', `N0', `Tpost', `N', `Tobs')
        mata: st_local("se", strofreal(se))
    }
    
    *--------------------------------------------------------------------------*
    *Display results and save results
    *--------------------------------------------------------------------------*
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
    tempvar ii m adoption trt
    egen `ii' = group(`2')
    qui xtset `ii' `3'
    local Tmin = r(tmin)         //t min
    local T    = r(tmax)         //t max
    local balanced = r(balanced) //balance

    if ("`balanced'"!="strongly balanced")==1 {
        di as err "Not Strongly Balanced data"
        exit 198
    }
	
    *Define treated periods
    bys `2' `4': egen `m' = min(`3')
    by  `2': egen `adoption' = max(`m')
    qui by `2': replace `adoption' = 0 if `adoption'==`Tmin'
    qui levelsof `adoption' if `adoption'>0, local(trt)
    qui tab `adoption' if `adoption'>0
    local length = r(r) 
    mata: Results = J(`length', 3, .)

    *filter data and SDiD
    mata: i = 1
    tempfile base
    qui save `base'

    ****
    keep `2' `ii' `adoption'
    bys `2': gen N=_n
    qui keep if N==1
    drop N
    rename `2' state
    rename `ii' statenumber
    rename `adoption' adoption
    sort adoption
    bys adoption: gen st=_n
    qui count if adoption==0
    local r = r(N)

    qui levelsof adoption if adoption>0, local(adop)
    foreach a of local adop {
        qui count if adoption==`a'
        local adop_`a'=r(N)
    }

    foreach a of local adop {
        local rr=`r'+1
        local ir=1
        while (`ir'<=`adop_`a'') {
            qui replace st=`rr' if adoption==`a' & st==`ir'
            local ++rr
            local ++ir
        }
    }
    tempfile resamplebase
    qui save `resamplebase'
	****

    use `base', clear
    foreach t of local trt {
        tempvar id id2 diff tr post_treat
        qui keep if `adoption'==`t' | `adoption'==0
        scalar time=`t'
        gen `post_treat' = 0
        qui replace `post_treat' = 1 if `3'>=`adoption' & `adoption'!=0
        qui sum `post_treat' if `post_treat'==1
        scalar obs=r(N)
        *----------------------------------------------------------------------*
        *- Create some temporal variables and locals                          -*
        *----------------------------------------------------------------------*
        egen `id' = group(`2')
        qui xtset `id' `3'
        local N_`t' = r(imax)                   //number of units
        qui tab `id' if `4'==1
        local Ntr_`t' = r(r)                    //number of treated units
        local N0_`t'  = `N_`t'' - `Ntr_`t''     //number of control units
        qui tab `3' if `3'>=`t'
        local Tpost_`t'  = r(r)                 //number of post times
        local T0_`t'     = `T'  - `Tpost_`t''   //max time of control
        local Tobs       = `T'  - `Tmin' +1     //number of times
        local Tpre_`t'   = `T0_`t'' - `Tmin' +1 //number of pre times
        local Ttrmin_`t' = `T0_`t'' + 1         //first year of treatment
        *----------------------------------------------------------------------*
        *- Calculate \zeta                                                    -*
        *----------------------------------------------------------------------*
        bys `id': egen `tr' = mean(`4')
        qui replace `tr' = 1 if `tr'!=0
        local EtaOmega  = (`Ntr_`t'' * `Tpost_`t'')^(1/4)
        local EtaLambda = 1e-6
        qui gen `diff' = `1' - L.`1'
        qui sum `diff' if `3'<=`T0_`t'' & `tr'==0
        local sig_`t' = r(sd)
        local ZetaOmega_`t'  = `EtaOmega'  * `sig_`t'' 
        local ZetaLambda_`t' = `EtaLambda' * `sig_`t''
        *----------------------------------------------------------------------*
        *- Preparing data                                                     -*
        *----------------------------------------------------------------------*
        tempfile data
        qui save "`data'"

        qui levelsof `3', local(times)                
        qui levelsof `3' if `3'<=`T0_`t'', local(timespre_`t') 
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
        drop Y1
        local i=2
        foreach n of local times {
            ren Y`i' t`n'
            local ++i
        }
        *All data for estimator
        mkmat _all, matrix(Yall_`t')
        mata: Yall_`t' = st_matrix("Yall_`t'")
        mata: Yall_`t' = Yall_`t'[1..`N_`t'',1..`Tobs']
        egen promt = rowmean(t`Ttrmin_`t''-t`T')
        drop t`Ttrmin_`t''-t`T'
        local r=`N_`t''+1
        qui set obs `r'

        forvalues tm=`Tmin'/`T0_`t'' {
            qui sum     t`tm' if id>`N0_`t''
            qui replace t`tm' = r(mean) in `r'
        }    

        qui drop if id>`N0_`t'' & id!=.
        mkmat _all, matrix(Y)
        *----------------------------------------------------------------------*
        *- Matrices for optimization                                          -*
        *----------------------------------------------------------------------*
        *Matrix A & b : Lambda 
        clear
        qui svmat Y, names(col)
    
        local vr `timespre_`t'' promt
        foreach tm of local vr {
            if "`tm'"=="promt" local n ""
            if "`tm'"!="promt" local n "t"
            qui sum `n'`tm' if id<=`N0_`t''
            qui replace `n'`tm' = `n'`tm' - r(mean) 
        }

        qui keep in 1/`N0_`t''
        mkmat promt, matrix(b_l)
        mkmat t`Tmin'-t`T0_`t'', matrix(A_l)
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

        local max=`N0_`t''+1
        forvalues tm=1/`max'  {
            qui sum t`tm'
            qui replace t`tm' = t`tm' - r(mean)
        }
			
        mkmat t`max', matrix(b_o)
        mkmat t1-t`N0_`t'', matrix(A_o)

        local mindec_`t'=(1e-5 * `sig_`t'')^2
        local col_l = colsof(A_l)
        local col_o = colsof(A_o)
        mat def lambda_l = J(1, `col_l', 1 / `col_l')
        mat def lambda_o = J(1, `col_o', 1 / `col_o')

        mata: A_l=st_matrix("A_l")
        mata: b_l=st_matrix("b_l")
        mata: A_o=st_matrix("A_o")
        mata: b_o=st_matrix("b_o")
        *--------------------------------------------------------------------------*
        *TAU
        *--------------------------------------------------------------------------*
        #delimit ;
        mata: tau = estTau(A_l,b_l,A_o,b_o,st_matrix("lambda_l"),st_matrix("lambda_o"),
                           Yall_`t', `ZetaLambda_`t'', `ZetaOmega_`t'',`mindec_`t'',
                           `Ntr_`t'',`N_`t'',`Tobs',`Tpost_`t'');
        #delimit cr

        ereturn matrix lambda_`t' lambda
        ereturn matrix omega_`t'  omega
        matrix LAMBDA_L_`t' = e(lambda_`t')
        matrix LAMBDA_O_`t' = e(omega_`t')
        mata: st_local("tau", strofreal(tau))

        mata: lambda_o_`t' = st_matrix("LAMBDA_O_`t'")
        mata: lambda_l_`t' = st_matrix("LAMBDA_L_`t'")

        mata: Results[i,1] = st_numscalar("time")
        mata: Results[i,2] = st_numscalar("obs")
        mata: Results[i,3] = tau
        mata: i = i+1
        use `base', clear
    }
    scalar drop time obs
    mata: Results[.,2] = Results[.,2]/sum(Results[.,2])
    mata: ATT = sum(Results[.,2] :* Results[.,3])
    mata: st_local("ATT", strofreal(ATT)) 	

    *--------------------------------------------------------------------------*
    *Standard error: bootstrap
    *--------------------------------------------------------------------------*
    if "`vce'"=="bootstrap" {
        set seed `seed'
        local b = 1
        local B = `reps'
        mata: ATT_b = J(`B', 1, .)
    
        /*if (`Ntr'==1)==1 {
            di as err "It is not possible to do Bootstrap se because there is only one treated unit"
            exit 198
        }*/
		
        dis "Bootstrap replications (`reps'). This may take some time."
        dis "----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5"

        clear
        while `b'<=`B' {
            use `resamplebase', clear
            bsample , cluster(statenumber) idcluster(bootState)
            tempfile resamplebase_b
            qui save `resamplebase_b'
            qui count if adoption == 0
            local r1 = r(N)
            qui count if adoption != 0
            local r2 = r(N)

            if (`r1'==0 | `r2'==0) {
                *di "all units are control or treated"
            }	
            else {
                display in smcl "." _continue
                if mod(`b',50)==0 dis "     `b'"
                qui levelsof adoption if adoption > 0, local(trt)
                qui tab adoption      if adoption > 0
                local length = r(r) 
                mata: results = J(`length', 6, .)

                mata: i = 1
                foreach t of local trt {
                    scalar time=`t'
                    qui keep if adoption==`t' | adoption==0
                    qui putmata ind1 = st, replace
                    qui putmata ind2 = st if adoption==0, replace

					mata: ind1 = sort(ind1,1)
					mata: ind2 = sort(ind2,1)
                    mata: Ntrb = length(ind1) - length(ind2) //treated units
                    mata: N = length(ind1)                   //N units
                    mata: st_local("Ntrb", strofreal(Ntrb)) 
                    mata: st_local("N", strofreal(N)) 	
                    mata: new_d=(Yall_`t'[ind1,],ind1)		
                    mata: st_matrix("new_d", new_d)

					clear
					qui svmat new_d, names(col)
                    local i=1
                    foreach n of local times {
                        ren c`i' t`n'
                        local ++i
                    }
                    local nv=`Tobs'+1
                    ren c`nv' id
                    sort id	
                    egen promt = rowmean(t`Ttrmin_`t''-t`T')
                    drop t`Ttrmin_`t''-t`T'					
                    local r=`N'+1
					
                    qui set obs `r'
                    forvalues tm=`Tmin'/`T0_`t'' {
                        qui sum     t`tm' if id>`N0_`t''
                        qui replace t`tm' = r(mean) in `r'
                    }    

                    qui drop if id>`N0_`t'' & id!=.
                    mkmat _all, matrix(Y)
                    *----------------------------------------------------------*
                    *- Matrices for optimization                              -*
                    *----------------------------------------------------------*
                    *Matrix A & b : Lambda 
                    clear
                    qui svmat Y, names(col)
                    local vr `timespre_`t'' promt
                    foreach tm of local vr {
                        if "`tm'"=="promt" local n ""
                        if "`tm'"!="promt" local n "t"
                        qui sum `n'`tm' if id<=`N0_`t''
                        qui replace `n'`tm' = `n'`tm' - r(mean) 
                    }

                    qui keep if id<`N0_`t''
                    mkmat promt, matrix(b_lb)
                    mkmat t`Tmin'-t`T0_`t'', matrix(A_lb)
                    local col_l = colsof(A_lb)
                    local row_l = rowsof(A_lb)
                    mata: A_lb = st_matrix("A_lb")
                    mata: b_lb = st_matrix("b_lb")
                    *Matrix A & b : Omega
                    clear
                    qui svmat Y, names(col)
                    drop promt id
                    gen id = _n
                    qui sum id
                    local max=`r(max)'
                    qui reshape long t, i(id) j(a)
                    qui reshape wide t, i(a) j(id)
                    drop a
                    local umax=`max'-1
                    forvalues tm=1/`umax'  {
                        qui sum t`tm'
                        qui replace t`tm' = t`tm' - r(mean)
                    }
			
                    mkmat t`max', matrix(b_ob)
                    mkmat t1-t`umax', matrix(A_ob)
                    mata: A_ob = st_matrix("A_ob")
                    mata: b_ob = st_matrix("b_ob")
					
                    mata: l_o = sum_norm(lambda_o_`t''[,ind2]) //actualizamos w omega
					mata: Yallb=Yall_`t'[ind1,]				

                    #delimit ;
                    mata: tau = estTau(A_lb, b_lb, A_ob, b_ob, 
                                        st_matrix("LAMBDA_L_`t'")',l_o, 
                                        Yallb, `ZetaLambda_`t'', `ZetaOmega_`t'',
                                        `mindec_`t'',Ntrb,`N_`t'',`Tobs',`Tpost_`t'');
                    #delimit cr
					
					scalar obs=`Ntrb'*`Tpost_`t''
					
                    mata: results[i,1] = st_numscalar("time")
                    mata: results[i,2] = st_numscalar("obs")
                    mata: results[i,3] = tau
                    mata: i = i+1
                    use `resamplebase_b', clear
                }
                scalar drop time obs
                mata: results[.,2] = results[.,2]/sum(results[.,2])
                mata: ATT_b[`b',] = sum(results[.,2] :* results[.,3])
                local ++b
            }
        }
        mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(ATT_b)))
        mata: st_local("se", strofreal(se))
    }
		
    *--------------------------------------------------------------------------*
    *Standard error: placebo
    *--------------------------------------------------------------------------*
    else if "`vce'"=="placebo" {
        dis "Standard error estimation under construction for staggered adoption"
    }
    *--------------------------------------------------------------------------*
    *Standard error: jackknife
    *--------------------------------------------------------------------------*
    else if "`vce'"=="jackknife" {
        dis "Standard error estimation under construction for staggered adoption"
    }
	
	*Display results
    di as text " "
    di as text "{c TLC}{hline 8}{c TT}{hline 11}{c TRC}"
    di as text "{c |} {bf: ATT}   {c |} " as result %9.5f `ATT'  as text " {c |}"
    di as text "{c |} {bf: se}    {c |} " as result %9.5f `se' as text " {c |}"
    di as text "{c BLC}{hline 8}{c BT}{hline 11}{c BRC}"   	
	
}

restore
end


*------------------------------------------------------------------------------*
*Mata functions
*------------------------------------------------------------------------------*
*Estimation
mata:
    real scalar estTau(matrix A_l, matrix b_l, matrix A_o, matrix b_o,
                       matrix lambda_l, matrix lambda_o, matrix Yall,
                       ZetaLambda, ZetaOmega, mindecrease, Ntr, N, Tobs, Tpost) {

        row_o = rows(A_o)
        row_l = rows(A_l)
        eta_o = row_o*ZetaOmega^2
        eta_l = row_l*ZetaLambda^2

        lambda_l = lambda(A_l,b_l,lambda_l,eta_l,ZetaLambda,100,mindecrease)
        lambda_l = sspar(lambda_l)
        lambda_l = lambda(A_l, b_l, lambda_l,eta_l,ZetaLambda, 10000,mindecrease)

        lambda_o = lambda(A_o, b_o, lambda_o,eta_o,ZetaOmega, 100,mindecrease)
        lambda_o = sspar(lambda_o)
        lambda_o = lambda(A_o, b_o, lambda_o,eta_o,ZetaOmega, 10000,mindecrease)

        st_matrix("lambda", lambda_l')
        st_matrix("omega",  lambda_o')


        tau = (-lambda_o, J(1,Ntr,1/Ntr))*Yall*(-lambda_l, J(1,Tpost,1/Tpost))'
        return(tau)
}
end
        
*minimization
mata:
real vector lambda(matrix A, matrix b, matrix x, eta, zeta, maxIter, mindecrease) {    
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
    real matrix sspar(matrix V)
{
    W = mm_cond(V :<= max(V)/4, 0, V)
    W = W :/ sum(W)
    return(W)
}
end
	
*sum normalize
mata:
    real vector sum_norm(matrix O)
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
    real matrix smerge(matrix A, matrix B)
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
    real vector jackknife(matrix Y, matrix L, matrix O, c, tp, N, T) {
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


			
