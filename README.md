This repository contains code and data to reproduce all figures from our paper:

Leonardo Novelli, Fatihcan M. Atay, Jürgen Jost, Joseph T. Lizier (2020).
[Deriving pairwise transfer entropy from network structure and motifs](https://arxiv.org/abs/1911.02931).
arXiv:1911.02931.

(More instructions in the `code` folder).

## Figures

### Figure 1
Generated manually to illustrate the motifs.


### Figure 2: Pairwise transfer entropy as a function of source and target in-degree.
![Results](figureOutputs/Fig2.png)
Figure 2. The pairwise transfer entropy (TE) increases with the in-degree of the source and decreases with the in-degree of the target, regardless of the sign of the link weights. The TE is plotted as a function of the source and target indegree. The results were obtained from 10 000 simulations of scale-free networks of 100 nodes generated via preferential attachment and the TE was averaged over all the node pairs with the same source and target in-degree. Note that the values in the lower left corner are the result of an average over many samples, since most of the node pairs have low in-degree. There are progressively fewer samples for higher in-degree pairs, and none for most pairs in the upper-right corner (absence indicated by the white colour).


### Figure 3: Average transfer entropy as a function of the rewiring probability in Watts-Strogatz ring networks.
![Results](figureOutputs/Fig3.png)
Figure 3. Average transfer entropy (TE) as a function of the rewiring probability in Watts-Strogatz ring networks. For positive link weights, the pairwise TE is higher in clustered networks than in random networks, due to the higher number of clustered motifs. For each value of the rewiring probability, the results for 10 simulations on different networks are presented (low-opacity markers) in addition to the mean values (solid markers). The plot shows that the approximation based on all the motifs up to order 4 (green curve) is closer to the theoretical values (orange curve) than the approximation based on the in-degrees and directed motifs alone (red curve) or on the directed motifs alone (violet curve). The empirical values are also shown (blue curve) as a validation of the theoretical results.


### Figure 4: Pairwise transfer entropy as a function of source and target in-degree in random Boolean networks.
![Results](figureOutputs/Fig4.png)
Figure 4. Pairwise transfer entropy (TE) as a function of the source and target in-degrees in random Boolean networks.
Similarly to the linear Gaussian case (Figure 2), the TE increases with the in-degree of the source and decreases with the in-degree of the target. The results were obtained from 10 000 simulations of scale-free networks of 100 nodes generated via preferential attachment. The TE was averaged over all the node pairs with the same in-degrees. The values in the lower left corner are the result of an average over many samples, since most of the node pairs have low in-degrees. There are progressively fewer observations for higher in-degrees and none in the upper-right corner (absence indicated by white colour).


### Figure 5: Average transfer entropy as a function of the rewiring probability in Watts-Strogatz ring networks with a random Boolean dynamics.
![Results](figureOutputs/Fig5.png)
Figure 5. Average empirical transfer entropy as a function of the rewiring probability in Watts-Strogatz ring networks with a random Boolean dynamics. The results for 10 simulations on different networks are presented (low-opacity markers) in addition to the mean values (solid markers).