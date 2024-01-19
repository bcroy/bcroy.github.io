---
layout: post
title:  "Eigendecompositions of random matrices for network analysis"
date:   2023-10-26 12:33:00 -0400
categories: RDPGs
---
<script src="/assets/js/mathjax-config.js" defer></script>

_You may have to refresh for the equation labels to display properly_

_This is a note for Prof. Alan Edelman. I have been working on this and related stuff for awhile, inspired by the work from [Carey Priebe](https://www.ams.jhu.edu/~priebe/) and colleagues at Johns Hopkins. In particular, I have been trying to come up with an alternative, more "geometric" approach to the paper by Mele et. al.[^Mele]. I think there's some novel stuff here that might be of interest to the network science community._

_The thing I am trying to figure out now is in the last section on "Embedding covariances..."_

* [this text is a placeholder, will be replaced with TOC by command below]
{:toc}

## Part 1 -- the simple case
Consider an Erdos-Renyi random graph $$G(n,p)$$, with $$n$$ nodes and the (symmetric) probability of an edge between them drawn according to a Bernoulli with probability $$p$$. Then, the adjacency matrix representation is an $$n \times n$$ symmetric matrix $$A$$, with $$A_{i,j} \in \{0,1\}$$ according to Bernoulli($$p$$).

Our goal here is to obtain a different representation of $$G$$ using the "random dot product graph" formalism. We seek embedding vectors $$x_i \in \mathbb{R}^{d}, d \leq n$$ for the nodes $$i$$ in $$G$$, such that $$A_{ij} = \vec{x_i} \cdot \vec{x_j}$$

_(I'm a fan of the RDP formulation -- it's simple and powerful and leads to some cool geometric interpretations as I hope to show.)_

We'll proceed with an example, generating an ER random graph with $$p=.16$$ and $$n=1000$$, and accomplish the embedding as follows.

Compute the eigendecomposition of $$A$$, to obtain $$A = Q \Lambda Q^T$$, with $$\lambda_i$$ ordered in decreasing order (and the columns of $$Q$$ ordered accordingly.)

Choose lower embedding dimension $$d$$, e.g. $$d=2$$ (we choose 2 since later we'll need 2 dimensions, though really we know there is only one parameter $$p$$ in this model.) Choosing $$d$$ is a fraught topic, though a common approach is to look for the "knee" in the eigenvalues, which we can see after the first eigenvalue. 

![Eigenvalues of ER(1000,.16)](/assets/images/ER_eigvals.png)

If we histogram the eigenvalues, we clearly see a single eigenvalue outside of the bulk.
In fact, the top 5 eigenvalues (in magnitude) are as follows:
{% highlight python %}
[161.29670666, -23.16983909, 23.11578288, -23.05396533, 22.71975887] 
{% endhighlight %}
We also have the very nice semi-circle law on display.
_I think we covered Wigner's semi-circle law in class... but unfortunately this is before I started coming in early October._ This paper by Guionnet[^guionnet] has more to say about Bernoulli random matrices.

![Eigenvalue histogram of ER(1000,.16)](/assets/images/ER_eigvals_hist.png)

Take the first $$d$$ columns of $$Q$$, yielding the $$n \times d$$ matrix $$Q_d$$, and the $$d \times d$$ diagonal matrix $$\Lambda_d$$ (i.e $$\lambda_1 \dots \lambda_d$$ on the diagonal.) Then, $$A_d = Q_d \Lambda_d Q_d^T$$.

Then let $$X = Q_d \Lambda_d^{\frac{1}{2}}$$. So, $$A_d = X X^T$$. For the moment, let's assume the largest $$d$$ eigenvalues are non-negative so the square root doesn't give any trouble. The extension to a "generalized RDPG" does allow for negative eigenvalues, see Rubin-Delanchy et. al. (2022)[^Rubin-Delanchy].

So, we'll use the set of $$n$$ $$d$$-dimensional row vectors $$X_i$$ as the "embeddings" of the nodes in $$G$$. With $$d=2$$, we visualize the node embeddings as seen in the scatter plot below.

![Embedding of nodes for ER(1000,.16) in d=2](/assets/images/ER_p=.16_scatter.png)

Wow, that's cool! What a nice, apparently axis-aligned ellipse. 
<!-- I colored the points $$x_i$$ by the number of non-zero entries in the correponding row. This corresponds to the "degree" of the node -- i.e. the number of incident edges. -->
What can we say about the position of this point cloud, and its variance? 

### Mean

Consider first the position. Imagine the "unreduced" embeddding (denote as $$X_n$$ for the full $$n$$ dimensions. Then $$X_ n = Q \Lambda^{\frac{1}{2}}$$, and $$\vec{x}, \vec{y}$$ as individual rows from $$X_n$$. We know that $$E\left[\vec{x} \cdot \vec{y}\right] = p$$, because $$\vec{x} \cdot \vec{y}$$ yields the corresponding entry in $$A$$, and the entries of $$A$$ have expected value $$p$$ since they are drawn from a Bernoulli with parameter $$p$$.

So, $$E\left[\vec{x} \cdot \vec{y}\right] = p$$,\\
$$E\left[\sum_{i=1}^n x_i y_i\right] = p$$,\\
$$\sum_{i=1}^n E [x_i y_i] = p$$,\\
$$\sum_{i=1}^n E [x_i] E[y_i] = p$$, (justified by independence)\\
$$\sum_{i=1}^n \mu_i^2 = p$$, (i think justified since $$\vec{x}$$ and $$\vec{y}$$ are drawn from same set of nodes and have same coordinate expectations)\\
$$\|\vec{\mu}\| = \sqrt{p}$$

Now... I am just going to assert that we'll extend this whole length of $$\sqrt{p}$$ along the first coordinate (i.e. the first eigenvector), so $$\mu_1 = \sqrt{p}$$. Thus, all the other coordinates $$\mu_2 \dots \mu_n$$ are 0. There must be a better way to do this... but in any case, $$\vec{\mu} = \left[\sqrt{p}, 0, \dots , 0\right]$$

In fact, the mean vector in 2-dimensions is [.40062, .00025], close to what we would expect. Hey, not bad!

### Covariance

Now consider the covariance of this point cloud. We compute the covariance matrix $$\mathrm{Cov}\left[X\right] = \frac{1}{n} \left(X-\vec{\mu}\right)^T \left(X-\vec{\mu}\right) = \frac{1}{n} X^T X - \vec{\mu}^T \vec{\mu}$$.

The left side becomes $$\frac{1}{n}X^T X = \frac{1}{n}\Lambda_d^{\frac{1}{2}} Q_d^T Q_d \Lambda_d^{\frac{1}{2}} = \frac{1}{n}\Lambda_d$$ since the columns of $$Q$$ are orthonormal.

The right side $$\vec{\mu}^T \vec{\mu}$$ becomes a matrix (this is an outer product since $$\vec{\mu}$$ is a row vector), with $$\mu_1^2 = p$$ for the top-left (1,1) element, and zero everywhere else.

In other words, the variance along the first dimension is $$\frac{\lambda_1}{n} - p$$, and the variance along the remaining dimensions are the eigenvalues $$\frac{\lambda_2}{n} \dots \frac{\lambda_n}{n}$$.

We can check this by simply computing the empirical covariance of the embedding vectors obtaining 
$$\mathrm{Cov}\left[X\right] = \left[\begin{array}{rr}
7.9714\times 10^{-4} & -9.9795\times10^{-5} \\
-9.9795\times10^{-5} & 2.3139\times10^{-2}\\
\end{array}\right]$$

and comparing the diagonal entries against $$\frac{\lambda_1}{n} - p = 0.00079714$$, $$\frac{\lambda_2}{n} = 0.02313892$$ 

Wow, that's cool!

### Getting fancier

All this involved generating $$A$$ and computing the eigenvalues, but it seems like we should be able to say something just from knowing that $$A$$ is a symmetric Bernoulli matrix with parameter $$p$$. Indeed, Füredi & Komlós[^furedi] tell us something about what to expect for $$\lambda_1$$ and the next largest (in magnitude) eigenvalue. 

Specifically, for our example $$\lambda_1$$ can be approximated by a normal distribution with expectation $$(n-1) p + 1$$ and variance $$2 p (1-p)$$, which we calculate as $$\mathcal{N}\left(160.84, 0.2688\right)$$ (see eqn 5, page 234). 

Equation 6 (page 234) states that the next largest eigenvalue (by absolute value) is upper bounded by $$2 \sqrt{n p (1-p)} + \mathcal{O}\left(n^\frac{1}{3}\log n\right)$$, and they note that by the semi-circle law it can't be much smaller than $$2 \sqrt{n p (1-p)}$$. Calculating these values does give a rather large approximate range of $$\left[23.186,  23.186 + \mathcal{O}\left(69.078\right)\right]$$. But for my purposes, the theoretical quantities are actually pretty good fits for the top 5 eigenvalues (repeated below).
{% highlight python %}
[161.29670666, -23.16983909, 23.11578288, -23.05396533, 22.71975887] 
{% endhighlight %}

For reference, there are other, more recent papers than Füredi & Komlós[^furedi] that dive deeper into the eigenvalues of random matrices. Specifically, the following series of recent papers: Benaych-Georges et. al. (2019)[^Benaych19], Benaych-Georges et. al. (2020)[^Benaych20] and Alt et. al. (2021)[^Alt]. However, I haven't really understood these yet, nor have I gone through the Füredi paper in detail beyond the basic result.

### Summary - Part 1

To sum up, computing node embeddings for an ER random graph using the eigendecomposition is well understood from a random matrix theory perspective. Given parameters $$p$$ and $$n$$ for an ER random graph, we can infer the expected mean and covariance of the node embedding point cloud. Conversely, given the $$n \times n$$ symmetric adjacency matrix $$A$$ from an observed ER random graph $$G$$, we can compute the node embeddings and infer the ER parameter $$p$$.

## Part 2 -- the thing I am trying to do

Here's where I am trying to take all of this.

Consider an Erdos-Renyi random graph $$G(n,p)$$, augmented with a node-level binary feature $$Z$$. The probability of an edge between two nodes that have the same binary feature is $$p+\beta$$, and an edge between nodes that have different binary features remains $$p$$. For example, the binary feature could denote a node's political party, and the value of $$\beta$$ reflecting how much it increases the edge probability between nodes with the same party affiliation.

As above, we are interested in the random dot product formulation and computing node embeddings. But now, we also wish to examine how the introduction of $$\beta$$ alters the geometry of the embeddings. 

We begin by simulating 3 networks with $$n=1000$$, $$p=.16$$ and $$\beta \in \left\{0,.15,.5\right\}$$. We assign each node's binary feature randomly with probability $$.5$$ and examine the effect of $$\beta$$. Their adjacency matrices (after permuting to group nodes by $$Z$$) are shown below. As expected, for $$\beta = 0$$ there is no discernable difference in the number of edges between nodes with $$Z=0$$ and $$Z=1$$. But as $$\beta$$ increases the block structure is more evident. (Note that we are being careful to choose $$\beta$$ such that $$p+\beta \leq 1$$ so that it remains a valid probability.)

![Adjacency matrices](/assets/images/Mele_p=.16,beta=0,.15,.5_adjmx.png)

Computing the eigenvalues of the corresponding adjacency matrices, we see that for $$\beta > 0$$, a second eigenvalue emerges.

![Eigenvalues](/assets/images/Mele_p=.16,beta=0,.15,.5_eigvals.png)
![Eigenvalue histogram](/assets/images/Mele_p=.16,beta=0,.15,.5_eigvals_hist.png)

We now compute and visualize the node embeddings in the same manner as the simple case above. That is, we compute the $$n \times d$$ dimensional node embedding matrix $$X = Q_d \Lambda_d^{\frac{1}{2}}$$ based on the eigendecomposition of the adjacency matrix. We show the node embeddings on the same plot for all three networks below. Since $$Z$$ is an observable node feature, we color the node embeddings accordingly. As expected, as $$\beta$$ increases the embedding clusters are better separated.

![Empirical clusters](/assets/images/empirical_clusters_p=.16,beta=0,.15,.5.png){:width="55%"}

### Geometric view of $$\beta$$
In this section, we show how the offset $$\beta$$ in "probability space" can be viewed geometrically, corresponding to the increasing separation of the clusters above.

We begin by showing how the offset to edge probabilty $$\beta$$ can be represented as a shift vector in embedding space $$\vec{a}$$, as follows. Let $$\vec{y}$$ correspond to the embedding vector for Democrats, and $$\vec{y}+\vec{a}$$ as the embedding vector for Republicans. The edge probabilities between nodes with each possible covariate pairing is as follows:

$$\begin{alignat}{3}
\mathsf{RR:} \ & p + \beta & = & \ \langle \vec{y}+\vec{a}, \vec{y}+\vec{a}\rangle    \label{eqRR}\\
\mathsf{RD:} \ & p         & = & \ \langle \vec{y}+\vec{a}, \vec{y}\rangle      \label{eqRD}\\
\mathsf{DD:} \ & p + \beta & = & \ \langle \vec{y}, \vec{y}\rangle        \label{eqDD}\\
\end{alignat}
$$

Then, we can combine equations (\ref{eqRD}) and (\ref{eqDD}) to get $$\vec{y} \cdot \vec{y} + \vec{y}\cdot\vec{a} + \beta = \vec{y}\cdot\vec{y}$$, implying

$$\beta = -\vec{y}\cdot\vec{a} \label{eqBeta1}$$

We can also combine equations (\ref{eqRR}) and (\ref{eqDD}) to obtain $$\vec{y}\cdot\vec{y} + 2\vec{y}\cdot\vec{a} + \vec{a}\cdot\vec{a} = \vec{y}\cdot\vec{y}$$, which simplifies to $$2\vec{y}\cdot\vec{a} + \vec{a}\cdot\vec{a} = 0$$, implying that 

$$\beta = \frac{\|\vec{a}\|^2}{2} \label{eqBeta2}$$

We can then rewrite equation (\ref{eqDD}) as 

$$p =\|\vec{y}\|^2 - \beta \label{eqP1}$$

### Inference
The above derivation leads to a method for estimating $$p$$ and $$\beta$$ from an observed adjacency matrix with node covariates. First, compute the node embedding vectors (the row vectors $$\vec{x}_i$$ of $$X$$) as before via the eigendecomposition. Since the node covariate is observable, then the mean vector $$\vec{x}_D$$ is simply $$\vec{x}_D = \textrm{avg}\left(\left\{\vec{x}_i: Z_i = 0\right\}\right)$$, and $$\vec{x}_R$$ computed similarly for $$Z_i=1$$.

The above derivation showed that the covariate could be represented as an offset vector $$\vec{a}$$, thus $$\hat{\vec{a}} = \vec{x}_R - \vec{x}_D$$. Applying equation (\ref{eqBeta2}) to $$\hat{\vec{a}}$$, we obtain an estimated $$\hat{\beta}$$ for each network. (Alternatively, equation (\ref{eqBeta1}) can also be used to estimate $$\beta$$ and gives similar results.) Using $$\vec{x}_D$$ as the estimate for $$\vec{y}$$ in equation (\ref{eqP1}) and the corresponding estimates for $$\beta$$, we obtain estimates for $$p$$ for each network. The table below collects these results.

| $$p$$ |$$\hat{p}$$| $$\beta$$ | $$\hat{\beta}$$          |
|---    |---|---        |---                       |
| .16   | 0.15993     | 0         | $$2.2075\times 10^{-5}$$ |
| .16   | 0.16076     | .15       | 0.15072                  |
| .16   | 0.162       | .5        | 0.49992                  |
{:.custom-table}

### Embedding positions as $$\beta$$ varies
The parameter $$\beta$$ governs the change in edge probability for nodes with the same covariate. For $$\beta = 0$$ the node covariates have no effect on the edge probabilities, and the graph reduces to a simple Erdos-Renyi random graph. On the other hand, as $$\beta$$ increases, two distinct clusters form, their means separated by offset vector $$\vec{a}$$ described above.

This is a little hand-wavy... but let's continue with the 2-dimensional case, and assume that $$\vec{a}$$ is axis aligned along the second eigenvector (and orthogonal to the first eigenvector). Then using equation (\ref{eqBeta2}) we know $$\vec{a} = \left[0, \sqrt{2\beta}\right]$$. 

Using $$\vec{a}$$ as above, and denoting the components of $$\vec{y}$$ as $$\left[y_0,y_1\right]$$, equation (\ref{eqBeta1}) yields $$\beta = -y_1 \sqrt{2\beta}$$, which gives $$y_1 = -\sqrt{\beta/2}$$.

We can solve for $$y_0$$ using equation (\ref{eqP1}), expanding as $$p = y_0^2 + \frac{\beta}{2} - \beta$$ and obtaining $$y_0 = \sqrt{p + \frac{\beta}{2}}$$.

Thus, the expected cluster centers are at

$$\vec{y} = \left[\sqrt{p + \frac{\beta}{2}}, -\sqrt{\frac{\beta}{2}}\right], \quad \vec{y} + \vec{a} = \left[\sqrt{p + \frac{\beta}{2}}, \sqrt{\frac{\beta}{2}}\right] \notag$$

We now rewrite this as a single equation. First, write the $$x$$ and $$y$$ coordinates as functions of $$\beta$$, so 

$$x^2(\beta) = p + \frac{\beta}{2}, \quad y^2(\beta) = \frac{\beta}{2} \notag$$

(_Note that squaring the $$y$$-coordinate for both $$\vec{y}$$ and $$\vec{y}+\vec{a}$$ yields $$\beta/2$$_). 
Combining these equations yields 

$$\frac{x^2(\beta) - y^2(\beta)}{p} = 1$$

This is the equation for a [hyperbola](https://en.wikipedia.org/wiki/Hyperbola). Now, applying this equation and plugging in the known values for $$p$$ and $$\beta$$ for our three example networks, we compute the expected cluster centers. The figure below overlays these predicted cluster centers on the actual clusters. In addition, we plot the underyling hyperbola curve. Looks like a pretty good fit! 

![Empirical clusters](/assets/images/empirical+theoretical_clusters_p=.16,beta=0,.15,.5.png){:width="55%"}



### Embedding ___covariances___ as $$\beta$$ varies

But now, this is where I am currently. What I really want is to derive the covariance matrices for the clusters. Given the mean position vectors and the corresponding covariance matrices, I could then say something about the expected errors for downstream inference tasks, which would give a more complete picture of this network model.

Here are some ideas that I am trying out:

Recall that $$X$$ is an $$n \times d$$ matrix of embeddings. Let $$X_1$$ be the $$n_1 \times d$$ sub-matrix corresponding to the subset of embeddings for cluster 1, and $$X_2$$ the $$n_2 \times d$$ sub-matrix for cluster 2, where $$n_1 + n_2 = n$$. As a simplification for now, it might come in handy later to assume $$n_1 = n_2 = \frac{1}{2}n$$... ideally we don't want to rely on this assumption though.

We wish to know $$\mathrm{Cov}\left[X_1\right]$$. If we can figure that out, then we can understand the orientation and shapes of the ellipses. As above, 

$$
\mathrm{Cov}\left[X_1\right] = \frac{1}{n_1} \left(X_1-\vec{\mu_1}\right)^T \left(X_1-\vec{\mu_1}\right) = \frac{1}{n_1} X_1^T X_1 - \vec{\mu_1}^T \vec{\mu_1} \notag
$$

and analgously for $$X_2$$. Expanding the left hand side, we have
$$X_1^T X_1 = \Lambda^\frac{1}{2}Q_1^T Q_1 \Lambda^\frac{1}{2}$$, where $$Q_1$$ corresponds to the $$n_1$$ rows of $$Q$$ for cluster 1. Note that $$Q_1^T Q_1$$ is not the identity matrix as was the case for $$Q^T Q$$, though it is symmetric. We write $$Q_1^T Q_1 =
\left[\begin{array}{rr}a_1 & b_1 \\
b_1 & c_1\\
\end{array}\right]$$. Then, $$\frac{1}{n_1}X_1^T X_1 = 
\frac{1}{n_1}\left[\begin{array}{cc}\lambda_1 a_1 & \sqrt{\lambda_1 \lambda_2} b_1 \\
\sqrt{\lambda_1 \lambda_2}b_1 & \lambda_2 c_1\\
\end{array}\right]$$

Taking the right hand side $$\vec{\mu_1}^T \vec{\mu_1}$$ and using the position vector $$\vec{\mu_1} = \left[\sqrt{p+\frac{\beta}{2}}, -\sqrt{\frac{\beta}{2}} \right]$$ derived above, we have
the outer product $$\vec{\mu_1}^T \vec{\mu_1} = 
\left[\begin{array}{cc}
p+\frac{\beta}{2} & -\sqrt{\frac{\beta}{2}\left(p+\frac{\beta}{2}\right)}\\
 -\sqrt{\frac{\beta}{2}\left(p+\frac{\beta}{2}\right)} & \frac{\beta}{2}
\end{array}\right]
$$.

Note that if $$\beta = 0$$, then right hand side matrix becomes $$\left[\begin{array}{cc}p& 0 \\ 0 & 0\end{array}\right]$$

It's a bit messy, but putting these together we have 

$$
\mathrm{Cov}\left[X_1\right] =
\frac{1}{n_1}\left[\begin{array}{cc}\lambda_1 a_1 & \sqrt{\lambda_1 \lambda_2} b_1 \\
\sqrt{\lambda_1 \lambda_2}b_1 & \lambda_2 c_1\\
\end{array}\right] -
\left[\begin{array}{cc}
p+\frac{\beta}{2} & -\sqrt{\frac{\beta}{2}\left(p+\frac{\beta}{2}\right)}\\
 -\sqrt{\frac{\beta}{2}\left(p+\frac{\beta}{2}\right)} & \frac{\beta}{2}
\end{array}\right]
\notag
$$

So now, if we could just figure the values for $$a_1, b_1, c_1$$ we would have a pretty complete characterization of this graph embedding!

Well, recall that $$a_1, b_1, c_1$$ are the elements of $$Q_1^T Q_1$$, and correspondly $$a_2, b_2, c_2$$ are the elements of $$Q_2^T Q_2$$, where $$Q_i$$ are the rows of $$Q$$ for group $$i \in \left\{1,2\right\}$$. It's not too hard to show that $$a_1 + a_2 = 1$$, $$c_1 + c_2 = 1$$, $$b_1 + b_2 = 0$$ since we know $$Q^T Q = I$$.

We can develop the relationship between $$\mathrm{Cov}\left[X_1\right], \mathrm{Cov}\left[X_2\right]$$ and $$\mathrm{Cov}\left[X\right]$$. I haven't done it in detail yet, but some things cancel but there is also a cross-term in $$\mathrm{Cov}\left[X\right]$$ that needs to be accounted for.

What if we focus on the rows of $$X_1$$ and $$X_2$$? We could rewrite every row $$i$$ of $$X_1$$ as $$\vec{\mu_1} + \vec{\epsilon_{1,i}}$$. To start, just consider one of the dimensions (say the "x" axis). Then $$\epsilon_i$$ is a mean 0 random variable, and maybe we could figure out its variance (since it's related to the Bernoulli). So what if we do this (taking $$\vec{q_1}$$ as the first column of $$Q_1$$ -- let's just deal with group 1 for now):\\
$$a_1 = \vec{q_1}^T \vec{q_1} = \sum_{i=1}^{n_1} q_{1,i}^2$$\
Let $$\bar{q_1} = \mathrm{mean}(q_1)$$, so $$q_{1,i} = \bar{q_1} + \epsilon_i$$.
Then $$a_1 = \sum_{i=1}^{n_1} \left(\bar{q_1} + \epsilon_i\right)^2$$ \
$$a_1 = n_1 \bar{q_1}^2 + \sum_{i=1}^{n_1} \epsilon_i^2$$ (the other terms with $$\sum_i \bar{q_1} \epsilon_i$$ vanish since $$\epsilon$$ is mean 0)
Now, we know something about $$\bar{q}$$ from earlier derivations. Do we know anything about the variance of $$\epsilon$$? 

Remember, $$\mu_1$$ and $$\epsilon_{1,i}$$ come directly from $$X$$, and a row of $$X$$ comes from the corresponding row of the adjacency matrix $$A$$ which really is just our set of $$n$$ Bernoulli trials. So... the variance ought to make its way through in some form.

### Ok... here's an idea
Following this development for PCA:
https://towardsdatascience.com/principal-component-analysis-part-1-the-different-formulations-6508f63a5553

We want to know the "variance of the projected data". In the sense that we are projecting each row of $$A$$ onto an eigenvector. For now let's not worry about the two groups, and just consider a single group with no offset $$\beta$$. So, start with row $$i$$ $$A_{i,\cdot}$$,  call it $$\vec{a}_i$$, and project onto the first eigenvector $$\vec{q_1}$$. Then $$\vec{x}_i = \frac{1}{\sqrt{\lambda_1}} \vec{a}_i \vec{q_1}$$. ...

Eventually, we can get an expression for the variance of $$\vec{x}$$ in terms of the variance in the original space, which will look something like 
$$\mathrm{Var}\left\[X_{\cdot,1}\right\] = \vec{q_1}^T S \vec{q_1}$$ where $$S$$ is the covariance matrix in the original space ... i.e. the space of our iid Bernoulli random variables. $$S$$ will be a symmetric $$n \times  n$$ matrix of covariances, so on the diagonal we will have $$\hat{sigma}^2 \approx p(1-p)$$. The off diagonals will have expected value of 0 (i think). Then what is the expected value of $$\vec{q_1}^T S \vec{q_1}$$? If we multiply it out, for each column of $$S$$ we would have a sum with many components of $$\vec{q_1}$$ multiplied by a mean 0 value and a single component of $$\vec{q_1}$$ multiplied by the diagonal element which is $$p(1-p)$$. So then, the whole thing will end up being $$\vec{q_1}_i^2 \hat{\sigma}_i^2$$. But I think the expected value of $$\hat{\sigma}_i^2 = p(1-p)$$ times the length $$\|\vec{q_1}\|^2 = 1$$. Or something. Actually, $$E\left[S\right] = \mathrm{diag}\left(p(1-p)\right)$$, right??

This might be a better derivation for the single block case (with no covariate) from part 1 above, since it might be able to express the variance in terms of the Bernoulli variance. But unfortunately I'm not sure if this can be applied to the case with a covariate, because in that case there are two different variances, with the first variance multiplied by a first chunk of the eigenvector and the second variance is multiplied by the second chunk. 

### Another thought
Ok, here's an idea. In the adjacency matrix, a row in the top half generally looks like a row in the bottom half except that the left and right groups are swapped (at least for the case where the groups sizes are equal). But in terms of how many 1s are in the top left block and the bottom right block, they are about the same, and the same is true for the top right and bottom left. Thus, I think we can argue that the top half of the values in the eigenvector are statistically similar to the bottom half of the values. Because if you grab a row from the bottom, swap left and right halves, and multiply it by the eigenvector it ought to yield a value that looks like the top half of the eigenvector. So that maybe gives us a constraint on the values in the eigenvector... they should be pretty similar, perhaps having the same mean and variance?

Interestingly, the squared norm of each half of the eigenvector is about the same proportion as the number of nodes in each group...


## Summary
Up to this point, I have focused on the simple ER random graph and the ER random graph with a single binary covariate. However, I have started on some work to show how this geometric interpretation can be applied to the stochastic block model. This is basically a "part 3" to this story -- see another draft writeup [here](http://bcroy.github.io/rdpgs/2023/06/12/mele-rdp.html#2-block-and-beyond-sbm-as-an-rdpg).

Questions:
- Can we determine the eigenvalues of a block symmetric matrix? In this case, the overall adjacency matrix has a symmetric block structure, each of which is a simple Bernoulli random matrix. Füredi & Komlós[^furedi] tells us something about the eigenvalues of Bernoulli random matrices, but what about the concatenation of them? Or maybe we could set it up as the sum of block matrices, and then can we say something about its eigenvalues based on the eigenvalues of the constituents?



# References
[^guionnet]: Guionnet, A. (2023). [Bernoulli Random Matrices.](https://ems.press/books/standalone/262/5172), also on arXiv [here](https://arxiv.org/abs/2112.05506)

[^ttao]: [Topics in random matrix theory](), T. Tao - 2023

[^furedi]: Füredi, Z., Komlós, J. (1981). [The eigenvalues of random symmetric matrices.](https://link.springer.com/article/10.1007/BF02579329) Combinatorica 1, 233–241.

[^Alt]: Alt, J., Ducatez, R., & Knowles, A. (2021). [Extremal eigenvalues of critical Erdős–Rényi graphs.](https://doi.org/10.1214/20-AOP1483) The Annals of Probability, 49(3), 1347–1401. 

[^Benaych19]: Benaych-Georges, F., Bordenave, C., & Knowles, A. (2019). [Largest Eigenvalues of Sparse Inhomogeneous Erdős–Rényi Graphs.](https://www.jstor.org/stable/26730468) The Annals of Probability, 47(3), 1653–1676.

[^Benaych20]: Benaych-Georges, F., Bordenave, C., & Knowles, A. (2020). [Spectral radii of sparse random matrices.](https://doi.org/10.1214/19-AIHP1033) Annales de l’Institut Henri Poincaré, Probabilités et Statistiques, 56(3), 2141–2161. 

[^Rubin-Delanchy]: Rubin-Delanchy, P., Cape, J., Tang, M., & Priebe, C. E. (2022). A Statistical Interpretation of Spectral Embedding: The Generalised Random Dot Product Graph. Journal of the Royal Statistical Society Series B: Statistical Methodology, 84(4), 1446–1473. https://doi.org/10.1111/rssb.12509

[^Mele]: Mele, A., Hao, L., Cape, J., & Priebe, C. E. (2021). Spectral inference for large Stochastic Blockmodels with nodal covariates (arXiv:1908.06438). arXiv. https://doi.org/10.48550/arXiv.1908.06438


<!--
Another way to do it..

in the text, use this [Guionnet](#guionnet)

to link within the page to the anchor below

[Bernoulli Random Matrices](https://ems.press/books/standalone/262/5172), also on arXiv [here](https://arxiv.org/abs/2112.05506), by Alice Guionnet (2023)
{: #guionnet}
-->

<!--Maybe relevant 
Distribution of Eigenvalues of Random Real Symmetric Block Matrices , https://arxiv.org/abs/1908.03834


-->
