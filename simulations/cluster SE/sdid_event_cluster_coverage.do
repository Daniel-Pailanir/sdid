clear

cap mata: mata drop cluster_data()
mata:
void cluster_data(q, n_q, T, alpha) {
    N = q*n_q
    M = J(N, 2 + T, .) // Cluster, Group, (Y_{g,t})_{t=1}^T
    M[1..N, 2] = (1..N)'
    for (j = 1; j <= q; j++) {
        M[((j-1)*n_q+1)..(j*n_q), 1] = J(n_q, 1, j)
        V = I(n_q) + alpha * uniform(1,1) :* (J(n_q,n_q,1) - I(n_q))
        M[((j-1)*n_q+1)..(j*n_q), 3..(T+2)] = matpowersym(V, 0.5)* rnormal(n_q,T,0,1)
    }
    st_store(.,.,M)
}
end

set seed 0
local CC = 100
local GC = 10
local TT = 10
local B = 100

mat res_sim = J(`B',28,.)
forv h = 1/`B' {
    qui {
        drop _all
        set obs `=`CC'*`GC''
        gen C = .
        gen G = .
        forv t = 1/`TT' {
            gen Y`t' = .
        }
        mata: cluster_data(`CC', `GC', `TT', 1)
        reshape long Y, i(C G) j(T)

        // No treatment effect, Basic Design, 2 treated clusters, 5 treated groups each
        gen D = C >= 90 & mod(G,2) == 0 & T>=`TT'/2
    }

    sdid_event Y G T D
    mat res_sim[`h',1] = (0>=e(H)[1,3]&0<=e(H)[1,4])
    forv j = 1/6 {
        mat res_sim[`h',1+`j'] = (0>=e(H)[1+`j',3]&0<=e(H)[1+`j',4])
    }

    sdid_event Y G T D, cluster(C)
    mat res_sim[`h',8] = (0>=e(H)[1,3]&0<=e(H)[1,4])
    forv j = 1/6 {
        mat res_sim[`h',8+`j'] = (0>=e(H)[1+`j',3]&0<=e(H)[1+`j',4])
    }

    sdid_event Y G T D, vce(placebo)
    mat res_sim[`h',15] = (0>=e(H)[1,3]&0<=e(H)[1,4])
    forv j = 1/6 {
        mat res_sim[`h',15+`j'] = (0>=e(H)[1+`j',3]&0<=e(H)[1+`j',4])
    }

    sdid_event Y G T D, vce(placebo) cluster(C)
    mat res_sim[`h',22] = (0>=e(H)[1,3]&0<=e(H)[1,4])
    forv j = 1/6 {
        mat res_sim[`h',22+`j'] = (0>=e(H)[1+`j',3]&0<=e(H)[1+`j',4])
    }

    di as err `h'
}

svmat res_sim
di "vce(bootstrap):"
forv j = 1/7 {
    qui sum res_sim`j'
    if `j' == 1 di "ATT: " r(mean) _continue
    else di "Effect `=`j'-1': " r(mean) _continue
    qui sum res_sim`=`j'+7'
    di " " r(mean)
}
di "vce(placebo):"
forv j = 1/7 {
    qui sum res_sim`=14+`j''
    if `j' == 1 di "ATT: " r(mean) _continue
    else di "Effect `=`j'-1': " r(mean) _continue
    qui sum res_sim`=`j'+21'
    di " " r(mean)
}