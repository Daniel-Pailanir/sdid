### Results from sdid_cluster_coverage.do

Coverage | Baseline | Cluster |
:-------- | :--------: | :--------: |
vce(bootstrap)  | 0.75  | 0.92  |
vce(placebo)  | 0.75  | 0.96  |
vce(jackknife)  | 0.75  | 0.95  |

_Simulation Details_: 100 simulations with periods $t = 1,..., T = 10$, clusters $c = 1, ..., C = 100$, each containing $N_c = 10$ groups each ($G = 1000$). Let $Y_{g,t}$ be the outcome of group $g$ at time $t$, $Y_g$ be the $1 \times T$ vector stacking the outcomes of group $g$, let $Y_c$ the $N_c \times T$ matrix stacking the vectors $Y_g$ for all $g$ in cluster $c$, and let $Y_{c,t}$ be the $t$-th column of $Y_c$ for $t = 1,..., T$. Let $\widetilde{Y}_{g,t} \sim N(0,1)$ be an i.i.d. collection for all $(g,t) \in \{1,..,G\} \times \{1,..,T\}$, and define $\widetilde{Y}_{g}$, $\widetilde{Y}_{c}$, $\widetilde{Y}_{c,t}$ as above. We set
$$
Y_c = V_c^{1/2} \widetilde{Y}_c
$$
where $V_c = I_{N_c} + U_c \left(1_{N_c}1'_{N_c} - I_{N_c}\right)$ is a $N_c \times N_c$ matrix with entries equal to 1 on the main diagonal, and equal to $U_c \sim U(0,1)$ elsewhere. This implies that 
$$
Y_{c,t} \sim N\left(\mathbf{0}_{N_c, 1}, V\right)
$$
i.e., within each time period, the outcomes of groups in cluster $c$ are a.s. positively correlated. 
As for the treatment rollout, we let 
$$
D_{g,t} = 1\{c_g \geq 90, g \bmod 2 = 0, t \geq 5\}
$$
Taking stock, we consider a non-staggered setting with 5 treated groups belonging to the same cluster.


