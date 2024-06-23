cap program drop sdid_event
program define sdid_event, eclass
syntax varlist(max = 4 min = 4 numeric) [if] [in] [, effects(integer 0) disag]

preserve
qui {
    marksample touse
    keep if `touse'
    keep `varlist'

    gen Y_XX = `1'
    gen G_XX = `2'
    gen T_XX = `3'
    gen D_XX = `4'

    sort G_XX T_XX 
    bys G_XX: egen ever_treated = max(D_XX)
    gen F_g_temp = D_XX*T_XX
    bys G_XX: egen C_temp = min(F_g_temp) if D_XX != 0
    replace C_temp = 0 if ever_treated == 0
    bys G_XX: egen C_XX = mean(C_temp)
    drop F_g_temp C_temp

    sum T_XX
    scalar T_max = r(max)
    gen T_g_XX = T_max - C_XX + 1 if ever_treated == 1
    sum T_g_XX
    local L_g = r(max)

    if `effects' == 0 local effects = `L_g'
    if `effects' > `L_g' {
        di as err "You have requested `effects' effects, but sdid_event can compute at most `L_g' dynamic effects."
        di as err "The estimation will resume with `L_g' effects."
        scalar effects = `L_g'
    } 
    else scalar effects = `effects'

    sdid Y_XX G_XX T_XX D_XX, vce(bootstrap) reps(2) mattitles
    mata: mata set matastrict off
    mata: omega = st_matrix("e(omega)")
    mata: lambda =st_matrix(" e(lambda)")
    matrix design = e(adoption)
    mata: st_numscalar("cohorts", rows(st_matrix("design")))

    mat res = J(1 + `L_g', `=cohorts', .)
    local rown "ATT_c"
    forv t = 1/`L_g' {
        local rown "`rown' Effect_c`t'"
    }
    local coln ""
    matrix c_weight = J(`=cohorts', 1, .)
    forv j = 1/`=cohorts' {
        local c = design[`j', 1]
        local coln "`coln' c=`c'"

        sum C_XX if C_XX == `c'
        scalar N_Post_`c' = T_max - r(mean) + 1
        unique G_XX if C_XX == `c'
        scalar N_Tr_`c' = r(unique)
        unique G_XX if inlist(C_XX, 0, `c')
        scalar N_`c' = r(unique)
        mat c_weight[`j', 1] = N_Tr_`c'

        mata: omega_temp = omega[1..(rows(omega) -1), `j']
        mata: lambda_temp = lambda[1..(rows(lambda) -1), `j']

        sort C_XX G_XX T_XX

        gen Y_XX_`c' = Y_XX if inlist(C_XX, 0, `c')
        mata: ATT_compute(st_data(., "Y_XX_`c'"), omega_temp, lambda_temp, st_numscalar("N_`c'"), st_numscalar("N_Post_`c'"), st_numscalar("N_Tr_`c'"))
        mat res[1, `j'] = ATT

        forv l = 1/`=N_Post_`c'' {
            gen Y_XX_`c'_`l' = Y_XX_`c' if (T_XX == `c' - 1 + `l' | T_XX < `c')
            mata: ATT_compute(st_data(., "Y_XX_`c'_`l'"), omega_temp, lambda_temp, st_numscalar("N_`c'"), 1, st_numscalar("N_Tr_`c'"))
            mat res[1 + `l', `j'] = ATT
        }
    }
    
    mata: ATT_aggte(st_matrix("res"), st_matrix("c_weight"), st_numscalar("effects"))

}   
    local rown_effects "ATT"
    forv j = 1/`=effects' {
        local rown_effects "`rown_effects' Effect_`j'"
    }
    mat rown H = `rown_effects'
    mat coln H = Estimate SE Switchers
    matlist H

    if "`disag'" != "" {
        mat res = res \ c_weight'
        local rown "`rown' Switchers"
        mat rown res = `rown' 
        mat coln res = `coln'
        di ""
        di "Disaggregated ATTs - Cohort level"
        matlist res
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
    ATT
    st_numscalar("ATT", ATT)
}
end

cap mata: mata drop ATT_aggte()
mata:
void ATT_aggte(B, W, E) {
    X = W' \ B
    H = J(rows(B), 3, .)
    for (j = 2; j <= rows(X); j++) {
        temp_mat = X[(1, j), .]
        temp_mat = select(temp_mat, temp_mat[2, .] :~= .)
        H[j-1, 1] = (temp_mat[1, .] :/ sum(temp_mat[1, .])) * temp_mat[2, .]'
        H[j-1, 3] = sum(temp_mat[1, .])
    }
    st_matrix("H",  H[1..(1+E), .])    
}
end