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
    reps(integer 50)
    controls(varlist numeric)
    graph(string)
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



*------------------------------------------------------------------------------*
*STAGGERED ADOPTION
*------------------------------------------------------------------------------*
**FOLLOWING 10 lines are the new implementation (so just need to bootstrap this)
**We could potentially move generation of treated and tyear into mata function...
tokenize `varlist'
tempvar treated ty tyear n
qui gen `ty' = `3' if `4'==1
qui bys `2': egen `treated' = max(`4')
qui by  `2': egen `tyear'   = min(`ty')
sort `3' `treated' `2'
gen `n'=_n
qui sum `3'
local mint=r(min)
qui putmata ori_id=`2' ori_pos=`n' if `3'==`mint' & `tyear'==., replace

if "`controls'"!="" unab conts: `controls'
**IDs (`2') go in twice here because we are not resampling
mata: data = st_data(.,("`1' `2' `2' `3' `4' `treated' `tyear' `conts'"))
mata: ATT = synthdid(data, 0, ., .)
mata: st_local("ATT", strofreal(ATT.Tau))
mata: OMEGA = ATT.Omega
mata: LAMBDA = ATT.Lambda


*some local
qui count if `3'==`mint'                 //total units
local N=r(N)
qui count if `treated'==0 & `3'==`mint' //control units
local co=r(N)
qui count if `treated'==1 & `3'==`mint' //treated units (total)
local tr=r(N)
local newtr=`co'-`tr'+1                 //start of treated units
qui tab `3'                             //T
local T=r(r)
mkmat `tyear' if `tyear'!=. & `3'==`mint', matrix(tryears) //save adoption values

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

    while `b'<=`B' {
        preserve
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
            mata: ATTB = synthdid(data,1,OMEGAB,LAMBDA)
            mata: ATT_b[`b',] = ATTB.Tau
            local ++b
        }
        restore
    }
    mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(ATT_b)))
    mata: st_local("se", strofreal(se)) 
}
    	
*--------------------------------------------------------------------------*
*Standard error: placebo
*--------------------------------------------------------------------------*
else if "`vce'"=="placebo" {
    set seed `seed'
    local b = 1
    local B = `reps'
    mata: ATT_p = J(`B', 1, .)
	
    /*if (`tr'>=`co')==1 {
        di as err "It is not possible to do Placebo se. Must have more controls than treated units."
        exit 198
    }*/
	
    dis "Placebo replications (`reps'). This may take some time."
    dis "----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5"
	
    while `b'<=`B' {
        preserve
        keep `1' `2' `3' `4' `tyear' `conts'
        tempvar r rand id
        qui gen `r'=runiform() in 1/`co'
        bys `2': egen `rand'=sum(`r')
        qui drop if `tyear'!=.        //drop treated units
        egen `id' = group(`rand')     //gen numeric order by runiform variable
		
        local i=1
        forvalues y=`newtr'/`co' {
            qui replace `tyear'=tryears[`i',1] if `id'==`y'
            local ++i
        }

		qui replace `4'=1 if `3'>=`tyear'
        bys `2': egen `treated' = max(`4')
		
        display in smcl "." _continue
        if mod(`b',50)==0 dis "     `b'"
			
        qui putmata psam_id=`2' if `tyear'==. & `3'==`mint', replace
        mata: indicator=smerge(psam_id, (ori_id, ori_pos))
        mata: OMEGAP=OMEGA[(indicator\rows(OMEGA)),]		
		
        mata: data = st_data(.,("`1' `2' `2' `3' `4' `treated' `tyear' `conts'"))
        mata: ATTP = synthdid(data,1,OMEGAP,LAMBDA)
        mata: ATT_p[`b',] = ATTP.Tau
		local ++b
        restore
    }
    mata: se = sqrt((`B'-1)/`B') * sqrt(variance(vec(ATT_p)))
    mata: st_local("se", strofreal(se))
}

*--------------------------------------------------------------------------*
*Standard error: jackknife
*--------------------------------------------------------------------------*
else if "`vce'"=="jackknife" {
    dis "Standard error estimation under construction for staggered adoption"
}
ereturn local se `se' 
ereturn local ATT `ATT'
    
*Display results
di as text " "
di as text "{c TLC}{hline 8}{c TT}{hline 11}{c TRC}"
di as text "{c |} {bf: ATT}   {c |} " as result %9.5f `ATT'  as text " {c |}"
di as text "{c |} {bf: se}    {c |} " as result %9.5f `se' as text " {c |}"
di as text "{c BLC}{hline 8}{c BT}{hline 11}{c BRC}"   	

end


*------------------------------------------------------------------------------*
*Mata functions
*------------------------------------------------------------------------------*
mata:
struct results {
    real matrix Omega
    real matrix Lambda
    real scalar Tau
}
end

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
    struct results scalar synthdid(matrix data, inference, OMEGA, LAMBDA) {
        data  = sort(data, (6,2,4))
        units = panelsetup(data,2)
        NT = panelstats(units)[,3]
        treat=panelsum(data[.,(2,5)],units)
        treat[,1]=treat[,1]/NT
        treat[,2]=1*(treat[,2]:>=1) + 0*(treat[,2]:==1)
        Nco = sum(data[,7]:==.)/NT
        controlID = uniqrows(select(data[.,2],data[,7]:==.))
        //controls = selectindex(treat[,2]:==1)
        
        
        //Adjust for controls in xysnth way
        if (cols(data)>7) {
            K = cols(data)
            cdat = data[selectindex(data[,6]:==0),(1,2,4,8..K)]
            cdat = select(cdat, rowmissing(cdat):==0)
            //CAN THIS BE OPTIMIZED FOR FIXED EFFECTS???
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
            // Estimate regression
            XX = quadcross(X,1 , X,1)
            Xy = quadcross(X,1 , y,0)
            beta = invsym(XX)*Xy
            beta = beta[1::NX]
            X = data[.,8::K]
            data[,1]=data[.,1]-X*beta
        }

        //Find years which change treatment
        trt = select(uniqrows(data[,7]),uniqrows(data[,7]):!=.)
        //Iterate over years, calculating each estimate
        tau    = J(rows(trt),1,.)
        tau_wt = J(1,rows(trt),.)

        if (inference==0) {
            Omega = J(Nco, rows(trt),.)
            Lambda = J(NT, rows(trt),.)
        }
        for(yr=1;yr<=rows(trt);++yr) {
            //trt[yr]
            cond1 = data[,7]:==trt[yr]
            cond2 = data[,7]:==.
            cond = cond1+cond2
            yNtr = sum(cond1)
            yNco = sum(cond2)
            //cond = cond + data[,6]:==.

            ydata = select(data,cond)
            yunits = panelsetup(ydata,2)
            yNT = panelstats(yunits)[,3]
            yNG = panelstats(yunits)[,1]
            yN  = panelstats(yunits)[,2]
            yNtr = yNtr/yNT
            yNco = yNco/yNT
            yTpost = max(ydata[,4])-trt[yr]+1

            pretreat = select(uniqrows(data[,4]),uniqrows(data[,4]):<trt[yr])
            Npre  = rows(pretreat)
            Npost = yNT-Npre

            //Calculate Zeta
            ndiff = yNG*(yNT-1)
            ylag = ydata[1..(yN-1),1]
            ylag = (. \ ylag)
            //ALL BELOW IS DOING THIS

            diff = ydata[.,1]-ylag

            first = mod(0::(yN-1),yNT):==0
            postt = ydata[,4]:>=trt[yr]
            dropc = first+postt+ydata[,6]
            prediff = select(diff,dropc:==0)
            sig_t = sqrt(variance(prediff))

            EtaLambda = 1e-6
            EtaOmega = (yNtr*yTpost)^(1/4)
            yZetaOmega  = EtaOmega*sig_t
            yZetaLambda = EtaLambda*sig_t
            //Generate Y matrices
            ids = uniqrows(ydata[.,2])
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
            A_o = Y0[,1..Npre]':-mean(Y0[,1..Npre]')
            //b_o = Y1[.,1..Npre]':-mean(Y1[.,1..Npre]')
            b_o = mean(Y1[.,1..Npre])':-mean(mean(Y1[.,1..Npre])')


            //Calculate Tau for t
            if (inference==1) {
                lambda_l = select(LAMBDA'[.,1::Npre],LAMBDA'[,rows(LAMBDA)]:==trt[yr])
                l_o = select(OMEGA'[.,1::yNco],OMEGA'[,rows(OMEGA)]:==trt[yr])
                lambda_o = sum_norm(l_o) //update so prior weights sum to 1
            }
            else {
                lambda_l = J(1,cols(A_l),1/cols(A_l))
                lambda_o = J(1,cols(A_o),1/cols(A_o))
            }
            mindec = (1e-5*sig_t)^2

            //Find optimal weight matrices
            eta_o = rows(A_o)*yZetaOmega^2
            eta_l = rows(A_l)*yZetaLambda^2
            lambda_l = lambda(A_l,b_l,lambda_l,eta_l,yZetaLambda,100,mindec)
            lambda_l = sspar(lambda_l)
            lambda_l = lambda(A_l, b_l, lambda_l,eta_l,yZetaLambda, 10000,mindec)

            lambda_o = lambda(A_o, b_o, lambda_o,eta_o,yZetaOmega, 100,mindec)
            lambda_o = sspar(lambda_o)
            lambda_o = lambda(A_o, b_o, lambda_o,eta_o,yZetaOmega, 10000,mindec)
            if (inference==0) {
                Lambda[.,yr] =  (lambda_l' \ J(Npost,1,.))
                Omega[.,yr] = lambda_o'
            }
            tau[yr] = (-lambda_o, J(1,yNtr,1/yNtr))*Y*(-lambda_l, J(1,Npost,1/Npost))'
            tau_wt[yr] = yNtr*Npost
        }
        if (inference==0) {
            Omega =  (Omega, controlID)
            Omega =  (Omega \ (trt', .))
            Lambda = (Lambda, uniqrows(data[,4]))
            Lambda = (Lambda \ (trt', .))
        }
        tau_wt = tau_wt/sum(tau_wt')
        ATT = tau_wt*tau
        struct results scalar r
        r.Omega = Omega
        r.Lambda = Lambda
        r.Tau = ATT
        return(r)
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

***CAN WE RE-WRITE mm_cond EASILY?  THIS WAY WE COULD REDUCE DEPENDENCIES...
*spar function
mata:
    real matrix sspar(matrix V) {
        W = mm_cond(V :<= max(V)/4, 0, V)
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

