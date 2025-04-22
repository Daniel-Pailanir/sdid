cap program drop sdid_event
program define sdid_event, eclass
syntax varlist(max = 4 min = 4) [if] [in] [, effects(integer 0) placebo(string) disag vce(string) brep(integer 50) method(string) covariates(string) vcov sb boot_ci combine(string)]
    version 12.0
    tempvar touse
    mark `touse' `if' `in'

    qui {

        foreach p in sdid unique {
            cap which `p'
            if _rc != 0 {
                noi ssc install `p', replace
            }
        }

        if "`method'" == "" local method "sdid"
        if !inlist("`method'", "sdid", "did", "sc") {
            di as err "Invalid syntax in method() option."
            exit
        }

        local m_sdid "Synthetic Difference-in-differences"
        local m_did "Difference-in-differences"
        local m_sc "Synthetic Control"

        tokenize `varlist'
        sort `2' `3'
        bys `2': egen ever_treated_XX = max(`4')

        if "`covariates'" != "" {
            gen Y_res_XX = `1'
            egen G_XX = group(`2')
            egen T_XX = group(`3')
            // Residualize Y on covariates, T and G FE on the sample of never-treated
            reg Y_res_XX `covariates' i.G_XX i.T_XX if ever_treated_XX == 0
            mat reg_res = e(b)
            local sel = 0
            foreach v in `covariates' {
                local sel = `sel' + 1
                replace Y_res_XX = Y_res_XX - reg_res[1, `sel'] * `v'
            }
            local varlist = subinstr("`varlist'", "`1'", "Y_res_XX", 1)
            drop G_XX T_XX
        }


        sdid_event_core `varlist' if `touse', effects(`effects') placebo(`placebo') method(`method') `disag' combine("`combine'")
        mat res_main = res
        mat H_main = H

        if "`vce'" == "" local vce "bootstrap"
        if !inlist("`vce'", "off", "bootstrap", "placebo") {
            di as err "Syntax error in vce option."
            di as err "Only {cmd:off}, {cmd:placebo} and {cmd:bootstrap} (dafalt) allowed."
            exit
        }
        if "`vce'" == "placebo" {
            unique `2' if ever_treated_XX == 0
            local N_co = r(unique)
            unique `2' if ever_treated_XX == 1
            if `N_co' < r(unique) {
                di as err "vce(placebo) can only be specified when the number of treated units is lower than the number of control units."
                exit
            }
        }

        if "`combine'" != "" {
            scalar tot_combine = length("`combine'") - length(subinstr("`combine'", ";", "", .)) + 1
            local arg_combine = "`combine';"

            mat H_cb = J(scalar(tot_combine), 3, .)
            local rown_b = ""
            forv k = 1/`=tot_combine' {
                local v_combine_`k' = substr("`arg_combine'", 1, strpos("`arg_combine'", ";")-1)
                local arg_combine = substr("`arg_combine'", strpos("`arg_combine'", ";") + 1, .)

                local t_k = 0
                foreach s in `v_combine_`k'' {
                    local t_k = `t_k' + 1
                }

                mat cset_`k' = J(`t_k', 3, .)
                local t_k = 1
                foreach s in `v_combine_`k'' {
                    mat cset_`k'[`t_k', 1] = `s'
                    mat cset_`k'[`t_k', 2] = H[`s' + 1, 1]
                    mat cset_`k'[`t_k', 3] = H[`s' + 1, 3]
                    local t_k = `t_k' + 1
                }

                mata: cbm = st_matrix("cset_`k'")
                mata: st_numscalar("cbn", sum(cbm[., 3]))
                mata: st_numscalar("cbt", (cbm[., 2]' * cbm[, 3])/sum(cbm[., 3]))

                mat H_cb[`k', 1] = scalar(cbt)
                mat H_cb[`k', 3] = scalar(cbn)
                local rown_b = "`rown_b' Cmb_Effect_`k'"
            }
            mat rown H_cb = `rown_b'
            mat coln H_cb = "Estimate" "SE" "Switchers"
        }

    }
    di as text "`m_`method''"

    if "`vce'" == "off" {
        matlist H_main, under

        local rownm: rown H_main
        mata: b = (st_matrix("H_main")[.,1])'
        mata: st_matrix("b", b)
        mat coln b = `rownm'

        ereturn post b
        ereturn matrix H = H_main
    }
    else {  
            qui {
                scalar breps = `brep'
                mata: boot_res_XX = J(st_numscalar("breps"), rows(st_matrix("H_main")), .)
                local failed = 0
            }

            di ""
            di "Boostrap replications (`brep'), `vce' mode."
            di "|0% " _dup(41) "-" " 100%|"
            di "|" _continue
            local counter = 1/50
            scalar r_XX = 1
            while r_XX <= breps {
                qui cap sdid_event_core `varlist' if `touse', effects(`effects') placebo(`placebo') method(`method') sampling(`vce') 
                if _rc == 0 {
                    mata: fail_coefs = rows(st_matrix("H_main")) - rows(st_matrix("H"))
                    mata: boot_res_XX[st_numscalar("r_XX"), .] = ((st_matrix("H")[., 1])', J(1, fail_coefs, .))
                    if (r_XX/`brep') > `counter' {
                        di "." _continue
                        local counter = `counter' + 1/50
                    }
                    scalar r_XX = r_XX + 1
                }
                else local failed = `failed' + 1
            }
            di "|" _newline

            {
                if "`boot_ci'" == "" {
                    mata: SE_add(boot_res_XX, st_matrix("H_main"))
                }
                else {
                    mata: SE_add_emp(boot_res_XX, st_matrix("H_main"))
                    di "CIs recovered from bootstrap distribution"
                }
            }

            local rownm: rown H_main
            mat rown H_SE = `rownm'
            mat coln H_SE = "Estimate" "SE" "LB CI" "UB CI" "Switchers"
            matlist H_SE, under

            mata: b = (st_matrix("H_SE")[.,1])'
            mata: V = diag((st_matrix("H_SE")[.,2]):^2)
            mata: st_matrix("b", b)
            mata: st_matrix("V", V)
            mat coln b = `rownm'
            mat rown V = `rownm'
            mat coln V = `rownm'

            ereturn post b V
            ereturn matrix H = H_SE

            if "`vcov'" != "" {
                mata: st_matrix("vcov", variance(boot_res_XX[.,2..(st_numscalar("effects")+1)]))
                local vcov_n = ""
                forv j = 1/`=effects' {
                    local vcov_n = "`vcov_n' Effect_`j'"
                }
                mat rown vcov = `vcov_n'
                mat coln vcov = `vcov_n'
                ereturn matrix vcov = vcov 
                if `=placebo' > 0 {
                    mata: st_matrix("vcov_pl", variance(boot_res_XX[.,(st_numscalar("effects")+2)..(st_numscalar("effects")+st_numscalar("placebo") + 1)]))
                    local vcov_n = ""
                    forv j = 1/`=placebo' {
                        local vcov_n = "`vcov_n' Placebo_`j'"
                    }
                    mat rown vcov_pl = `vcov_n'
                    mat coln vcov_pl = `vcov_n'
                    ereturn matrix vcov_pl = vcov_pl 
                }
            }

            if `failed' > 0 {
                di as result ""
                di as result "WARNING: Restarted `failed' bootstrap run(s) with no treated or control groups."
            }

            if "`combine'" != "" {
                mata: boot_res_cb_XX = J(rows(boot_res_XX), st_numscalar("tot_combine"), 0)
                forv k = 1/`=tot_combine' {
                    scalar cbid = `k'
                    scalar cbtw = H_cb[`k', 3]
                    local t_k = 1
                    foreach s in `v_combine_`k'' {
                        scalar id = `s'
                        scalar cbw = cset_`k'[`t_k', 3]
                        mata: boot_res_cb_XX[., st_numscalar("cbid")] = boot_res_cb_XX[., st_numscalar("cbid")] + boot_res_XX[., st_numscalar("id") + 1] * st_numscalar("cbw")
                        local t_k = `t_k' + 1
                    }
                    mata: boot_res_cb_XX[., st_numscalar("cbid")] = boot_res_cb_XX[., st_numscalar("cbid")] / st_numscalar("cbtw")
                }

                {
                    if "`boot_ci'" == "" {
                        mata: SE_add_cb(boot_res_cb_XX, st_matrix("H_cb"))
                    }
                    else {
                        mata: SE_add_cb_emp(boot_res_cb_XX, st_matrix("H_cb"))
                    }
                }
                
                mat H_cb = H_cb[.,1..1], H_cb_SE[.,2..5]
                mat coln H_cb = "Estimate" "SE" "LB CI" "UB CI" "Switchers"
            }
    }

    if "`combine'" != "" {
        di ""
        di "Combined ATTs"
        matlist H_cb, under
        ereturn matrix H_comb = H_cb
    }


    if "`disag'" != "" {
        di ""
        di "Disaggregated ATTs - Cohort level"
        matlist res_main, under
        ereturn matrix H_c = res_main
    }

    if "`sb'" != "" {
        mata: st_matrix("boot_res", boot_res_XX)
        ereturn matrix sb = boot_res
    }

    cap drop *_XX
end

cap program drop sdid_event_core
program define sdid_event_core, eclass
syntax varlist(max = 4 min = 4) [if] [in], effects(integer) method(string) [disag placebo(string) sampling(string) combine(string)]
preserve
qui {

    keep `if'
    sort `2' `3'
    gen Y_XX = `1'
    egen G_XX = group(`2')
    egen T_XX = group(`3')
    gen D_XX = `4'
    keep Y_XX D_XX G_XX T_XX ever_treated_XX

    if "`sampling'" != "" {
        if "`sampling'" == "bootstrap" {
            bsample, cluster(G_XX)
            sort G_XX T_XX
            bys G_XX T_XX: gen ID_XX = _n
            egen G_temp_XX = group(G_XX ID_XX)
            replace G_XX = G_temp_XX
            sort G_XX T_XX
            drop G_temp_XX ID_XX
        }
        if "`sampling'" == "placebo" {
            sort G_XX T_XX
            mata: treat = st_data(., ("G_XX", "T_XX", "D_XX", "ever_treated_XX"))
            mata: treat = select(treat, treat[., 4])[., 3]
            mata: st_matrix("treat", treat)
            keep if ever_treated_XX == 0

            bsample, cluster(G_XX)
            sort G_XX T_XX
            bys G_XX T_XX: gen ID_XX = _n
            egen G_temp_XX = group(G_XX ID_XX)
            replace G_XX = G_temp_XX
            sort G_XX T_XX
            drop G_temp_XX ID_XX

            drop D_XX ever_treated_XX
            gen id_r = uniform()
            bys G_XX: egen id_rg = mean(id_r)
            sort id_rg G_XX T_XX
            svmat treat
            rename treat1 D_XX
            replace D_XX = 0 if missing(D_XX)
            
            sort G_XX T_XX
            bys G_XX: egen ever_treated_XX = max(D_XX)            
        }
    }

    gen F_g_temp = D_XX*T_XX
    bys G_XX: egen C_temp = min(F_g_temp) if D_XX != 0
    replace C_temp = 0 if ever_treated_XX == 0
    bys G_XX: egen C_XX = mean(C_temp)
    drop F_g_temp C_temp

    sum T_XX
    scalar T_max = r(max)
    scalar T_min = r(min)

    // Number of feasible event study coefficients
    gen T_g_XX = T_max - C_XX + 1 if ever_treated_XX == 1
    sum T_g_XX
    local L_g = r(max)
    scalar L_g = `L_g'

    // Number of feasible placebo estimates
    sum C_XX if ever_treated_XX == 1
    local L_pl_g = r(max) - T_min

    if `effects' == 0 local effects = `L_g'
    if `effects' > `L_g' {
        di as err "You have requested `effects' effects, but sdid_event can compute at most `L_g' dynamic effects."
        di as err "The estimation will resume with `L_g' effects."
        scalar effects = `L_g'
    } 
    else scalar effects = `effects'

    scalar placebo = 0
    if "`placebo'" != "" {
        if "`placebo'" == "all" {
            scalar placebo = `L_pl_g'
        }
        else {
            if `placebo' > `L_pl_g' {
                di as err "You have requested `placebo' placebo, but sdid_event can compute at most `L_pl_g' placebo."
                di as err "The estimation will resume with `L_pl_g' placebo."
                scalar placebo = `L_pl_g'
            }
            else {
                scalar placebo = `placebo'
            }
        }
    }

    sdid Y_XX G_XX T_XX D_XX, vce(noinference) method(`method')
    mata: mata set matastrict off
    mata: omega = st_matrix("e(omega)")
    mata: lambda = st_matrix("e(lambda)")
    matrix design = e(adoption)
    mata: st_numscalar("cohorts", rows(st_matrix("design")))

    mat res = J(1 + `L_g', `=cohorts', .)
    local rown "ATT_c"
    forv t = 1/`L_g' {
        local rown "`rown' Effect_c`t'"
    }
    if "`placebo'" != "" {
        mat res_pl = J(`L_pl_g', `=cohorts', .)
        forv t = 1/`L_pl_g' {
            local rown "`rown' Placebo_c`t'"
        }
    }
    local coln ""
    matrix c_weight = J(`=cohorts', 1, .)
    matrix t_weight = J(`=cohorts', 1, .)
    forv j = 1/`=cohorts' {
        local c = design[`j', 1]
        local coln "`coln' c=`c'"

        sum C_XX if C_XX == `c'
        scalar N_Post_`c' = T_max - r(mean) + 1
        scalar N_Pre_`c' = r(mean) - T_min
        unique G_XX if C_XX == `c'
        scalar N_Tr_`c' = r(unique)
        unique G_XX if inlist(C_XX, 0, `c')
        scalar N_`c' = r(unique)
        mat c_weight[`j', 1] = N_Tr_`c'
        mat t_weight[`j', 1] = N_Post_`c' * N_Tr_`c'

        mata: omega_temp = omega[1..(rows(omega) -1), `j']
        mata: lambda_temp = lambda[1..(rows(lambda) -1), `j']

        sort C_XX G_XX T_XX

        gen Y_XX_`c' = Y_XX if inlist(C_XX, 0, `c')
        mata: mata set matastrict off
        mata: ATT_compute(st_data(., "Y_XX_`c'"), omega_temp, lambda_temp, st_numscalar("N_`c'"), st_numscalar("N_Post_`c'"), st_numscalar("N_Tr_`c'"))
        mat res[1, `j'] = ATT

        forv l = 1/`=N_Post_`c'' {
            gen Y_XX_`c'_`l' = Y_XX_`c' if (T_XX == `c' - 1 + `l' | T_XX < `c')
            mata: mata set matastrict off
            mata: ATT_compute(st_data(., "Y_XX_`c'_`l'"), omega_temp, lambda_temp, st_numscalar("N_`c'"), 1, st_numscalar("N_Tr_`c'"))
            mat res[1 + `l', `j'] = ATT
        }

        if "`placebo'" != "" {
            gen Y_XX_`c'_pre = Y_XX_`c' if (T_XX < `c')
            forv l = 1/`=N_Pre_`c'' {
                scalar pl = `l'
                mata: mata set matastrict off
                mata: ATT_compute_pl(st_data(., "Y_XX_`c'_pre"), omega_temp, lambda_temp, st_numscalar("N_`c'"), st_numscalar("pl"),st_numscalar("N_Tr_`c'"))
                mat res_pl[`l', `j'] = ATT
            }
        }
    }

    if "`placebo'" != "" {
        mat res = res \ res_pl
    }
    
    //scalar tot_est = scalar(effects) + scalar(placebo)
    mata: ATT_aggte(st_matrix("res"), st_matrix("t_weight"), st_numscalar("effects"), st_numscalar("placebo"), st_numscalar("L_g"), st_matrix("c_weight"))

}   

    local rown_effects "ATT"
    forv j = 1/`=effects' {
        local rown_effects "`rown_effects' Effect_`j'"
    }
    if "`placebo'" != "" {
        forv j = 1/`=placebo' {
            local rown_effects "`rown_effects' Placebo_`j'"
        }
    }
    mat rown H = `rown_effects'
    mat coln H = Estimate SE Switchers

    if "`disag'" != "" {
        mat res = res \ c_weight'
        local rown "`rown' Switchers"
        mat rown res = `rown' 
        mat coln res = `coln'
    }

restore
end

cap mata: mata drop ATT_compute()
mata:
void ATT_compute(Y, omega, lambda, G_max, N_Post, N_Tr) {
    Y_nm = select(Y, Y[.,1] :~= .)
    lambda_nm = select(lambda, lambda[.,1] :~= .)
    Y_mat = rowshape(Y_nm, G_max)
    ATT = ((-omega'), J(1, N_Tr, 1/N_Tr)) * Y_mat * ((-lambda_nm)\J(N_Post, 1, 1/N_Post))
    st_numscalar("ATT", ATT)
}
end

cap mata: mata drop ATT_compute_pl()
mata:
void ATT_compute_pl(Y, omega, lambda, G_max, pl, N_Tr) {
    Y_nm = select(Y, Y[.,1] :~= .)
    lambda_nm = select(lambda, lambda[.,1] :~= .)
    Y_mat = rowshape(Y_nm, G_max)
    Y_mat = Y_mat, Y_mat[.,cols(Y_mat)-pl+1]
    ATT = ((-omega'), J(1, N_Tr, 1/N_Tr)) * Y_mat * ((-lambda_nm)\1)
    st_numscalar("ATT", ATT)
}
end

cap mata: mata drop ATT_aggte()
mata:
void ATT_aggte(B, W, E, P, TE, S) {
    X = W' \ S' \ B
    H = J(rows(B), 3, .)
    for (j = 3; j <= rows(X); j++) {
        if (j == 3) {
            temp_mat = X[(1, j), .]
        }
        else {
            temp_mat = X[(2, j), .]
        }
        temp_mat = select(temp_mat, temp_mat[2, .] :~= .)
        H[j-2, 1] = (temp_mat[1, .] :/ sum(temp_mat[1, .])) * temp_mat[2, .]'
        H[j-2, 3] = sum(S[1..cols(temp_mat),.])
    }
    if (P == 0) {
        st_matrix("H",  H[1..(1+E), .])    
    }
    else {
        st_matrix("H",  H[(1..(1+E),(2+TE)..(1+TE+P)), .])    
    }
}
end

cap mata: mata drop SE_add()
mata:
void SE_add(V, B) {
    H = B[,1], J(rows(B), 3, .), B[, 3]
    for (j = 1; j <= rows(B); j++) {
        H[j, 2] = sqrt(variance(V)[j,j])
        H[j, 3] = H[j, 1] - 1.96 * H[j, 2]
        H[j, 4] = H[j, 1] + 1.96 * H[j, 2]
    }
    st_matrix("H_SE", H)
}
end

cap mata: mata drop SE_add_emp()
mata:
void SE_add_emp(V, B) {
    H = B[,1], J(rows(B), 3, .), B[, 3]
    for (j = 1; j <= rows(B); j++) {
        H[j, 2] = sqrt(variance(V)[j,j])
        V_nm = select(V[.,j], V[.,j]:~= .)
        K = sort(V_nm,1), ((1.. rows(V_nm))/rows(V_nm))'
        LT = sort(select(K, K[.,2]:<=0.025), 1) 
        if (rows(LT) == 0) {
            LT = K[1,.]
        }
        RT = sort(select(K, K[.,2]:>=0.975), 1)
        if (rows(RT) == 0) {
            RT = K[rows(K),.]
        }
        H[j, 3] = LT[rows(LT), 1]
        H[j, 4] = RT[1,1]
    }
    st_matrix("H_SE", H)
}
end

cap mata: mata drop SE_add_cb()
mata:
void SE_add_cb(V, B) {
    H = B[,1], J(rows(B), 3, .), B[, 3]
    for (j = 1; j <= rows(B); j++) {
        H[j, 2] = sqrt(variance(V)[j,j])
        H[j, 3] = H[j, 1] - 1.96 * H[j, 2]
        H[j, 4] = H[j, 1] + 1.96 * H[j, 2]
    }
    st_matrix("H_cb_SE", H)
}
end

cap mata: mata drop SE_add_cb_emp()
mata:
void SE_add_cb_emp(V, B) {
    H = B[,1], J(rows(B), 3, .), B[, 3]
    for (j = 1; j <= rows(B); j++) {
        H[j, 2] = sqrt(variance(V)[j,j])
        V_nm = select(V[.,j], V[.,j]:~= .)
        K = sort(V_nm,1), ((1.. rows(V_nm))/rows(V_nm))'
        LT = sort(select(K, K[.,2]:<=0.025), 1) 
        if (rows(LT) == 0) {
            LT = K[1,.]
        }
        RT = sort(select(K, K[.,2]:>=0.975), 1)
        if (rows(RT) == 0) {
            RT = K[rows(K),.]
        }
        H[j, 3] = LT[rows(LT), 1]
        H[j, 4] = RT[1,1]
    }
    st_matrix("H_cb_SE", H)
}
end