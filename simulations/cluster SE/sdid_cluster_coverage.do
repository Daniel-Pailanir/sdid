clear
qui do "sdid.ado"

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

mat res = J(`B',4,.)
forv h = 1/`B' {
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
    gen D = ((G>85&G<=90)|(G>95&G<=100)) & T>=`TT'/2
    sdid Y G T D, vce(bootstrap) 
    mat res[`h', 1] = (0>=e(b)[1,1]-1.96*sqrt(e(V)[1,1])&0<=e(b)[1,1]+1.96*sqrt(e(V)[1,1]))
    sdid Y G T D, vce(bootstrap) cluster(C)
    mat res[`h', 2] = (0>=e(b)[1,1]-1.96*sqrt(e(V)[1,1])&0<=e(b)[1,1]+1.96*sqrt(e(V)[1,1]))
    sdid Y G T D, vce(placebo) 
    mat res[`h', 3] = (0>=e(b)[1,1]-1.96*sqrt(e(V)[1,1])&0<=e(b)[1,1]+1.96*sqrt(e(V)[1,1]))
    sdid Y G T D, vce(placebo) cluster(C)
    mat res[`h', 4] = (0>=e(b)[1,1]-1.96*sqrt(e(V)[1,1])&0<=e(b)[1,1]+1.96*sqrt(e(V)[1,1]))
}

svmat res
forv j = 1/4 {
    qui sum res`j'
    di r(mean)
}