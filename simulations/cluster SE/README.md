### Results from sdid_cluster_coverage.do

Coverage | sdid Y G T D | sdid Y G T D, cluster(C) |
:-------- | :--------: | :--------: |
vce(bootstrap)  | 0.75  | 0.92  |
vce(placebo)  | 0.75  | 0.96  |
vce(jackknife)  | 0.75  | 0.95  |

### Simulation Details 
100 simulations with periods $t = 1,..., T = 10$, clusters $c = 1, ..., C = 100$, each containing $N\_c = 10$ groups each ($G = 1000$). Let $Y\_{g,t}$ be the outcome of group $g$ at time $t$, $Y\_g$ be the $1 \times T$ vector stacking the outcomes of group $g$, let $Y\_c$ the $N\_c \times T$ matrix stacking the vectors $Y\_g$ for all $g$ in cluster $c$, and let $Y\_{c,t}$ be the $t$-th column of $Y\_c$ for $t = 1,..., T$. 

Let $\widetilde{Y}\_{g,t} \sim N(0,1)$ be an i.i.d. collection for all $(g,t) \in \{1,..,G\} \times \{1,..,T\}$, and define $\widetilde{Y}\_{g}$, $\widetilde{Y}\_{c}$, $\widetilde{Y}\_{c,t}$ as above. We set 

$$Y\_c = V\_c^{1/2} \widetilde{Y}\_c$$

where $V\_c = I\_{N\_c} + U\_c \left(1\_{N\_c}1'\_{N\_c} - I\_{N\_c}\right)$ is a $N\_c \times N\_c$ matrix with entries equal to 1 on the main diagonal, and equal to $U\_c \sim U(0,1)$ elsewhere. This implies that 

$$Y\_{c,t} \sim N\left(\mathbf{0}\_{N\_c \times 1}, V\right)$$

i.e., within each time period, the outcomes of groups in cluster $c$ are a.s. positively correlated. 
As for the treatment rollout, we let 

$$D\_{g,t} = 1\lbrace c\_g \geq 90, g \bmod 2 = 0, t \geq 5 \rbrace $$

Taking stock, we consider a non-staggered setting with 5 treated groups belonging to the same cluster. Since the outcome is independent of the treatment, the true ATT is 0.


