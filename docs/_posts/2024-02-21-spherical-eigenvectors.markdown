---
layout: post
title:  "Spherical eigenvectors"
date:   2024-02-21 16:35:00 -0500
categories: Random matrices
---
<script src="/assets/js/mathjax-config.js" defer></script>

_You may have to refresh for the equation labels to display properly_

* [this text is a placeholder, will be replaced with TOC by command below]
{:toc}

## Introduction

I've been working on graph embeddings, focused on node embeddings for the Erdős-Renyi random graph and extending to embeddings for nodes that have binary covariates, inspired by the work of Mele et. al. (2021).[^Mele]

The purpose of this post is to collect some facts... and working assumptions... relevant to the stuff I've been doing with graph embeddings and random matrices.

## Random Dot Product Graphs

The Random Dot Product Graph (RDPG) views the adjacency matrix of a graph as arising from the dot products of underlying node embedding vectors. The method uses the eigendecomposition of the adjacency matrix to obtain vector representations for the nodes, such that the inner product of embedding vectors approximately reconstructs the original adjacency matrix.

To begin, we briefly introduce the formalism for an Erdős-Renyi (ER) random graph $$G(n,p)$$ characterized by $$n$$ nodes with probability $$p$$ of an edge between any pair of nodes. The corresponding $$n \times n$$  adjacency matrix $$A$$ is symmetric, with entries $$A_{ij} = A_{ji} \sim \mathrm{Bernoulli}(p)$$. 

Given an observed network (and observed adjacency matrix $$A$$), node embeddings can be defined as follows. Compute the eigendecomposition $$A = Q \Lambda Q^T$$, where all matrices are $$n \times n$$ and $$\Lambda$$ is a diagonal matrix consisting of eigenvalues $$\lambda_i$$ in decreasing order. By taking only the first $$d$$ columns of $$\Lambda$$ (corresponding to the columns containing the $$d$$ largest eigenvalues), one obtains the $$n \times d$$ dimensional node embeddings $$X = Q \Lambda_d^{1/2}$$. We are interested in the geometry of the $$d$$-dimensional node embeddings in $$X$$.

## Distribution of eigenvectors on the sphere
There are a number of interesting and relevant results that can help us in exploring the geometry of $$X$$. First, let's begin by examining the eigenvectors, the $$n$$ orthonormal columns of $$Q$$. Since each $$\|\vec{q}_i\| = 1$$ then $$\left\{\vec{q}_i\right\}$$ correspond to a set of $$n$$ points on $$\mathcal{S}^{n-1}$$. How are they distributed on the hypersphere? 

Goldstein (2017)[^Goldstein] tells us that such an orthonormal basis is expected to be uniformly distributed on the sphere. Put another way, any region $$\mathcal{R}$$ on the sphere, with probability measure $$u(\mathcal{R})$$, is expected to contain that fraction of the $$\vec{q}_i$$ points. A related result in Cai et. al. (2013)[^Cai] states that vectors corresponding to random points on the unit sphere are almost always nearly orthogonal. (I recall from Jim Anderson's neural networks class in college the rule of thumb that if you randomly generate a set of high dimensional unit vectors they'll be pretty close to orthogonal, so I guess this is where that comes from.) In any case, the Goldstein result is particularly useful because we are given orthgonal vectors (the $$n$$ eigenvectors of $$A$$), and it's useful to know that they should be distributed uniformly on the unit hypersphere. Going the other way, if we had a set of random $$n$$-dimensional vectors, which by Cai et. al. are expected to be nearly orthgonal and uniformly distributed, then presumably there is an orthonormal transformation that could (nearly) align them to the eigenvectors, right? By the way, see Liu, 2021[^Liu], section 3 for an informal but concise description of the relevant Goldstein theorem. In the _Experiments_ section below I perform some simulations to test the null hypothesis that the eigenvectors are uniformly distributed on the sphere.

## Expected norm of subvectors

Ok, so if we believe the eigenvectors of $$A$$ are uniformly distributed on the sphere, then a fact described in O'Rourke et. al. (2016)[^ORourke] puts it to use for our problem. Specifically, Theorem 2.3, pg. 6, regarding the distribution of mass. To put it in my own words, it says that for any $$n$$-dimensional vector $$\vec{v}$$ uniformly distributed on the unit sphere, computing the squared norm of only a subset $$S$$ of its coordinates yields a Beta distribution with mean 
$$\frac{|S|}{n}$$. This is written compactly as $$\|\vec{v}\|^2_S \sim \mathrm{Beta}\left(\frac{|S|}{2},\frac{n-|S|}{2}\right)$$. I think this is pretty interesting, and should be pretty easy to explore empirically.

But how would this be helpful for our purposes? Well, consider the scenario where the nodes are partitioned into two groups. I want to show that ...

__A concern:__ I am a little worried that lurking in the background of the O'Rourke paper, and papers by Van Vu and others, is the focus on Wigner matrices. I am pretty sure in the O'Rourke paper the expected norm over a subset of coordinates is generalized to any vector uniformly distributed on the sphere. Similarly, the Goldstein result seems to be general to any set of orthonormal vectors and makes no mention of Wigner matrices. But nonetheless, the adjacency matrices I am working on are not exactly Wigner matrices, and so need to be careful here.

## Experiments


### Uniform distribution on the sphere
It turns out there is a whole literature related to the statistics of distributions on the sphere. YARH! (**Y**et **A**nother **R**abbit **H**ole I went down...) There are also software packages, mostly in R, for testing distributions of points on the sphere $$S^{n-1}$$. The [`sphunif`](https://github.com/egarpor/sphunif) package seemed particularly good, along with the companion reference [^Garcia-Portugues]. Other references included

Some other references:
- Recent advances in directional statistics, by Arthur Pewsey and Eduardo García-Portugués
- CircStats R package: https://cran.r-project.org/web/packages/CircStats/CircStats.pdf
- Method to generate points: https://mathoverflow.net/questions/24688/efficiently-sampling-points-uniformly-from-the-surface-of-an-n-sphere
- 

#### Points that _should_ be uniformly distributed
We begin with a simple simulation -- generate a set of points that _should_ be uniformly distributed, and confirm that they are. One method to generate a set of points uniformly distributed on the sphere $$S^{n-1}$$ is, for each point, draw $$n$$ samples from a standard normal, collect into a vector, and normalize to unit length. For convenience, the `sphunif` provides the `r_unif_sph` function to generate uniformly distributed points on the sphere. The `unif_test` function provides multiple testing methods. I used the Rayleigh test and the CCF09 method, documented in [^Garcia-Portugues]. I generated 1000 vectors of dimension 1000 distributed uniformly on the sphere, testing with both the Rayleigh and CCF09[^Cuesta-Albertos] method, which did not reject the null hypothesis of uniformity. So that's good. Another interesting sanity check is to compute the pairwise dot products, to find that there are two regimes -- only 1000 dot products close to 1 (corresponding to the dot products of random vectors with themselves), and the remainder in a band $$[-.2,.2]$$ centered around 0. Specifically, for those in this band the mean and standard deviation are $$\hat{\mu} = 4.87005e^{-5}, \hat{\sigma}^2\approx.001$$. So, we see that a set of random vectors on the sphere $$S^{n-1}$$ are nearly orthogonal. Following the Stack Exchange post [here](https://stats.stackexchange.com/questions/85916/distribution-of-scalar-products-of-two-random-unit-vectors-in-d-dimensions)[^se1], the correlation is expected to be Beta distributed (but close to normal for large dimension $$n$$), with variance $$1/n$$, as observed.

#### Eigenvectors of the ER random graph

Next, take the eigenvectors of the ER random graph, which according to Goldstein et. al.[^Goldstein] should be uniformly distributed as well. Here I generated the ER random graph $$G(n=500,p=.16)$$ and computed the eigendecomposition of the adjacency matrix. I applied the same tests as above to the eigenvectors, once again failing to reject the null hypothesis that the eigenvectors are uniformly distributed on the sphere. Here of course, it makes less sense to histogram the dot products since they explicitly form an orthonormal basis and should either be 0 or 1. So far so good...

#### Eigenvectors of the ER random graph with binary covariate

We now consider the eigenvectors of the kinds of graphs we're interested in from the work in Mele et. al. (2021)[^Mele], which have different connection probabilities depending on the node covariates. In these graphs, each node is characterized by an observable binary covariate $$Z$$ assigned with probability $$q$$. The probability of an edge between a pair of nodes $$i,j$$ is $$p$$ when $$Z_i \neq Z_j$$, otherwise it is $$p + \beta$$. In other words, when nodes have the same covariate value (e.g. political party), the probability of forming an edge is increased by $$\beta$$. This yields a block structure in the adjacency matrix. So, for graphs generated with different choices of $$q$$ (assignment of covariate labels) and $$\beta$$ (effect of the covariate), how are the eigenvector distributions affected? It's less clear that these vectors should be uniformly distributed, although Goldstein (2017)[^Goldstein] suggests this should be the case.

I generated such graphs for the same base edge connection probability $$p=.16$$, and for several choices of $$q$$ and $$\beta$$. I repeated the same statistical tests on the eigenvectors assoicated with the graphs for each parameter setting. As before, the statistical tests do not seem to discern any non-uniformity to their distribution on the sphere.

Of course, in all cases involving eigenvectors of symmetric matrices, they are going to be mutually orthogonal and equidistant on the sphere, so they wouldn't exhibit any "clumping" that would be evidence for any non-uniformity. On the other hand, they aren't "uniform" in that they are constrained to be orthogonal to each other. So perhaps this isn't really an appropriate thing to test. 


But, the reason I am interested in this is due to the theorem described in O'Rourke in the section above on the expected norm of subvectors. That is, if the vectors are uniformly distributed, then the norm of a subset of their coordinates should behave in a predictable way. I am specifically interested in whether the norm of the subset of coordinates corresponding to nodes with a particular group behaves as described.

### Norms of subvectors

Lets now explore the result stated in O'Rourke (2016)[^ORourke], that a subvector consisting of $$S$$ coordinates of a vector uniformly distributed on the sphere has squared norm that is Beta distributed with mean 
$$\frac{|S|}{n}$$.

We begin with the approach above, generating vectors uniformly distributed on the sphere, sampling subsets of their coordinates, and computing the norms.

*Show the graph of the norms of random vectors as we vary 
$$|S|$$*

Next, we consider the eigenvectors. First we take random subsets of their coordinates as above, and note that the norms of those subvectors behave in the same way

*Show the graph of the norms of the eigenvectors as we vary 
$$|S|$$*

Finally, we consider the eigenvectors again, this time choosing the specific subset of coordinates associated with each group (i.e. $$Z=0$$ or $$Z=1$$)

But what about the distribution of angles of these subvectors... what are the dot products of $$\vec{e}_{1,S}$$ and $$\vec{e}_{2,S}$$? This is relevant because that's part of the cross-term in the covariance matrices.


### Dot products

The squared norm of the subvectors of eigenvectors, which we explored above, seem to follow the result from O'Rourke. How about the dot products between subvectors of different eigenvectors? Here we see something interesting... for almost all eigenvectors, the dot products of subsets (randomly chosen or partitioned by the covariate) are distributed around 0. Meaning, they are approximately orthogonal. However, not so for the "meaningful" eigenvectors we use for our embedding, corresponding to those with the largest eigenvalues. In fact, $$\vec{e}_{1,\mathcal{S}} \cdot \vec{e}_{2,\mathcal{S}}$$ is roughly $$-.5$$ for one subset and $$.5$$ for the complement. (They are constrained to sum to 0, since the full vectors are orthogonal.)


Some relevant posts on this topic:

[Distribution of scalar products of two random unit vectors in 𝐷 dimensions](https://stats.stackexchange.com/questions/85916/distribution-of-scalar-products-of-two-random-unit-vectors-in-d-dimensions)
Distribution of dot product of unit vectors uniformly distributed on the sphere... this has the answer, and the derivations. Main issue is that I'm not sure how applicable for subvectors...?

[Distribution of dot product of two unit random vectors](https://mathoverflow.net/questions/208937/distribution-of-dot-product-of-two-unit-random-vectors)
Less useful

[Distribution of Dot-Product of Two Independent Multivariate Gaussian Vectors](https://math.stackexchange.com/questions/3173291/distribution-of-dot-product-of-two-independent-multivariate-gaussian-vectors)
Distribution of dot products, but where vectors are not unit vectors. This has a good development by the question author... but the response is not especially tractable with stuff about Bessel functions, etc.

## References

[^ORourke]: O’Rourke, S., Vu, V., & Wang, K. (2016). Eigenvectors of random matrices: A survey. Journal of Combinatorial Theory, Series A, 144, 361–442. https://doi.org/10.1016/j.jcta.2016.06.008

[^Mele]: Mele, A., Hao, L., Cape, J., & Priebe, C. E. (2021). Spectral inference for large Stochastic Blockmodels with nodal covariates (arXiv:1908.06438). arXiv. https://doi.org/10.48550/arXiv.1908.06438

[^Furedi]: Füredi, Z., Komlós, J. (1981). [The eigenvalues of random symmetric matrices.](https://link.springer.com/article/10.1007/BF02579329) Combinatorica 1, 233–241.

<!-- Stuff on distribution of points on the sphere -->

[^Goldstein]: Goldstein, S., Lebowitz, J. L., Tumulka, R., & Zanghî, N. (2017). Any orthonormal basis in high dimension is uniformly distributed over the sphere. Annales de l’Institut Henri Poincaré, Probabilités et Statistiques, 53(2), 701–717. https://doi.org/10.1214/15-AIHP732

[^Cai]: Cai, T., Fan, J., & Jiang, T. (2013). Distributions of Angles in Random Packing on Spheres (arXiv:1306.0256). arXiv. https://doi.org/10.48550/arXiv.1306.0256

[^Liu]: Liu, W., Lin, R., Liu, Z., Xiong, L., Schölkopf, B., & Weller, A. (2021, March). Learning with hyperspherical uniformity. In International Conference On Artificial Intelligence and Statistics (pp. 1180-1188). PMLR.

[^Garcia-Portugues]: García-Portugués, E., & Verdebout, T. (2018). An overview of uniformity tests on the hypersphere (arXiv:1804.00286). arXiv. http://arxiv.org/abs/1804.00286

[^Cuesta-Albertos]: Cuesta-Albertos, J. A., Cuevas, A., and Fraiman, R. (2009). On projection-based tests for directional and compositional data. Stat. Comput., 19(4):367–380.

[^se1]: Distribution of scalar products of two random unit vectors in $$D$$ dimensions. https://stats.stackexchange.com/questions/85916/distribution-of-scalar-products-of-two-random-unit-vectors-in-d-dimensions
