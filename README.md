# A Bayesian CART for Modeling Gene Regulation Using Multiple High-Throughput Sequencing Data


#### Table of Contents
* [Motivation](#motivation)
* [Bayesian CART](#bayesian-cart)
* [Our Modification](#our-modification)
* [Implementation Details](#implementation-details)
* [Final Comments](#final-comments)


## Motivation

The question we are concerned with is as follows:
> If we have a bunch of RNASeq data of, say *E. coli*, and some ChIPSeq data of transcription factors (TFs) binding to genes’ promoter, how can we model the regulatory relationship between the TFs and the genes they bind?

There are some useful information and some caveats for us to consider:

1. RNASeq data reveals the **activity level** of TFs and genes in their mRNA form. This is useful, as we may expect some “correlation” between a TF and a gene when the TF is a true regulator of the gene. For example, if TF X is an activator for gene Y, we may expect a high X corresponding to a high Y.

2. Caveats for point 1: TF regulates gene in the form of protein, not mRNA, so it is really just an approximation in using RNASeq data to represent the TF’s activity level. For TFs that require post-translational modification to become active, this approximation will be bad.

3.	ChIPSeq data tells us which TF binds where at what affinity level. This is useful, as we may expect that when TF X is a true regulator of gene Y, X tends to bind somewhere near Y, with high affinity. 

4.	Caveats for point 2: not all binding are regulatory. This is easy to understand, if we think of binding as a thermodynamic event. But perhaps we need to be even more careful when we say that regulatory bindings are those with high affinity. There are some notable exceptions, but I need to look into my books to say anything concrete.


When we think about this problem, we imagine a modeling scheme that can roughly be divided into 2 steps:

* For a target gene, we want to sample the “true” TFs, i.e. regulating TFs as a subset of candidate TFs. The candidate TFs would be all TFs that bind near the gene. The sampling may have a weight proportional to the binding affinity, so that TFs with high binding affinity are more likely to be sampled. Moreover, in simple organism like *E. coli*, there shouldn’t be more than a couple of TFs being “true” TF for a gene, since we do not expect complex regulation in simple organisms.
* Once we have sampled the “true” TFs, we want to establish a “correlation” between the “true” TFs and the gene. In other words, we want something that we can conclude with, say, “when TF A is high and TF B is low, the gene is high; when TF A is low and TF B…”.

This sounds awfully like an MCMC sampling scheme coupled with a binary tree. But at the time we don’t know how to probabilistically represent a tree. So, we looked around and found a paper of about my age (I hope) that provided a good starting point for us:

> Hugh A. Chipman , Edward I. George & Robert E. McCulloch (1998) Bayesian CART Model Search, Journal of the American Statistical Association, 93:443, 935-948, DOI: 10.1080/01621459.1998.10473750

## Bayesian CART

The method tries to build a classification and regression tree (CART) in a Bayesian approach. More specifically, a binary tree is represented probabilistically. The data are grouped according to the binary decisions into leaves where the response variables are evaluated based on classification/regression task. Here is an example of such tree from the paper:

![Fig 2 of Chipman paper](/img/img1.png)

In our case, we can choose one TF at each decision node, and group samples into leaves. At each leaf, target genes’ data can fit to some parametric family of distributions, say Gaussian. 

The interesting part is how the model is specified. The prior for the tree is specified as:

![f](http://chart.apis.google.com/chart?cht=tx&chl=$P(\Theta,T)=P(\Theta|T)P(T)$)

Where T represents the tree and ![f](http://chart.apis.google.com/chart?cht=tx&chl=$\Theta$) the parameters at each leaf. On P(T), it is defined as follows: starting from a node ![f](http://chart.apis.google.com/chart?cht=tx&chl=$\eta$), if we are to grow it by creating two leaves, the probability of it splitting and assigning some splitting rule is:

![f](http://chart.apis.google.com/chart?cht=tx&chl=$P_{split}(\eta,T)P_{rule}(\rho|\eta,T)$)

where ![f](http://chart.apis.google.com/chart?cht=tx&chl=$P_{split}(\eta,T)=\alpha(1%2Bd_\eta)^\beta$). ![f](http://chart.apis.google.com/chart?cht=tx&chl=$d_\eta$) is the depth of node ![f](http://chart.apis.google.com/chart?cht=tx&chl=$\eta$), ![f](http://chart.apis.google.com/chart?cht=tx&chl=$\alpha$) and ![f](http://chart.apis.google.com/chart?cht=tx&chl=$\beta$) are some specified hyperparameters. ![f](http://chart.apis.google.com/chart?cht=tx&chl=$P_{rule}(\rho|\eta,T)$) can be specified by the user. We define it as a discrete distribution over all the candidate TFs (with weights being their binding affinity to the gene) and a discrete uniform distribution over all data values of the chosen TF.

On ![f](http://chart.apis.google.com/chart?cht=tx&chl=$P(\Theta|T)$), that is the parameters prior in each leaf, the authors proposed to use conjugate prior so that parameters can be marginalized in the joint distribution, and sampling can be simplified:

![f](http://chart.apis.google.com/chart?cht=tx&chl=$P(Y|X,T)={\int}P(Y|X,\Theta,T)P(\Theta|T)d\Theta$)

So in the paper when data in each leaf are assumed to have a normal distribution with unknown mean and prior, we can use a normal-inverse-gamma prior for the model parameters. The integration is pretty standard and can be found in many Bayesian statistics textbooks, such as *Bayesian Data Analysis* by Gelman *et al.*.

And the posterior of the tree is simply:

![f](http://chart.apis.google.com/chart?cht=tx&chl=$P(T|X,Y){\propto}P(Y|X,T)P(T)$)

To use MCMC to sample from the posterior, the authors proposed a few ways to propose new trees as the sampling kernel, but we have only implemented the following 3: 
* **Grow**: randomly pick a leaf node and split. In our case, we pick a leaf node and sample a TF to split it. If that TF has appeared in the ancestors of the leaf, we record down the splitting value at that ancestor, and sample a splitting value at the leaf different from that ancestral value. Also, we do not pick value that will result in a leaf with less than n data samples. Default setting of n is 20.
* **Prune**: close two leaf nodes and make their parent a new leaf.
* **Change**: pick an internal node and reassign a splitting rule. Similar to what we implement in Grow step, we ensure the resampled TF and splitting value have not occurred in the nodes among ancestors and offspring, and resampling will not result in leaf with too few samples.


## Our Modification

Clearly, we can specify the likelihood function at the leaves however we want, and if we use a conjugate prior for its parameters, we hardly need to change the whole scheme at all. A simple modification we did was to have a linear regression with Gaussian noise at each leaf. The regressors are all the ancestors of the leaf. Statistically there is little change, as we still specify a normal-inverse-gamma prior and derive the marginal likelihood in closed form. But now we have a very neat biological interpretation for the linear regression: since a leaf is a region specified by binary decisions of the tree, say “TF A is high and TF B is low”, the regression coefficient of TF A and TF B tells us their relationship with the target gene in that confined region.

Because of the outrageously poor support of LaTeX here, we don't show the derivation details. But again, there are many resources on deriving the marginal likelihood of Bayesian linear model. As an example, we found this derivation pretty clear:

> http://www.biostat.umn.edu/~ph7440/pubh7440/BayesianLinearModelGoryDetails.pdf


## Implementation Details

Some comments on our implementation. We created a `Tree` object which consists of 4 containers.Maps (hash table in Matlab). As an example here, T_max is a `Tree` object.

```matlab
>> T_max

T_max = 

  Tree with properties:

    lchild_ids: [7×1 containers.Map]
    rchild_ids: [7×1 containers.Map]
    parent_ids: [7×1 containers.Map]
         nodes: [7×1 containers.Map]
```

Each of them has node ID (integer) as keys. `parent_ids` has nodes’ parent’s node ID as values; `lchild_ids` has nodes’ left child node ID as values; `rchild_ids` has nodes’ right child node ID as values; `nodes` has `Node` objects as values. So if you want to know the ID of all leaf nodes, you can do:

```matlab
>> allIDs = cell2mat(T_max.lchild_ids.keys);
>> leafIDs = allIDs(isnan(cell2mat(T_max.lchild_ids.values)))

leafIDs =

     4     5     6     7
```

`Node` object stores information about each node. For example, the node ID, decision feature (i.e. TF ID), decision value (i.e. value for binary decision), ID of data that have passed through this node, etc. So if you want to know the row ID of your dataset (i.e. experiment sample ID) that made it through the leaf node 4, you can do:

```matlab
>> T_max.nodes(4).sample_ids;
```

The `main` function takes input the data file, which should include:

* nxp matrix of TF level, where n is number of data samples and p the number of candidate TFs.
* 1xp cell array of TF names.
* nxq matrix of gene level, where q is the number of equivalent genes (e.g. 3 genes that are in the same transcription unit and are regarded to be regulated exactly the same way, so you may want to use all 3 as response variable).
* 1xq cell array of gene names.
* 1xp vector of weights of TFs, the prior for the binding affinity of TFs.

We have a `data.mat` file in the repo as an example input.

Also in `main` function, there are some variables you can have fun with. `BLR` is a boolean that can set to true if to use Bayesian linear regression at each leaf; false if to use a Gaussian model described in Chipman’s paper. `n` is the smallest number of data samples in any leaf node, `alpha` and `beta` the hyperparameters for ![f](http://chart.apis.google.com/chart?cht=tx&chl=$P_{split}(\eta,T)$), and `iter` the max number of MCMC iterations. You can read more in the inline comments of the code.

By the way, while you are having fun with the code, if the process is interrupted by say an error, you may want to type `clear all` in the command to clear all objects before you start a new run, even if in the Matlab workspace window there is no objects. The issue may have to do with containers.Map being a reference type object, and when creating a new tree, if the old tree's objects are not cleared, the constructor of new tree will crumble.   

`main` function calls `MCMC` function and save the *tree with highest posterior probability* along with other outputs to a file in the format `beta_<beta hyperparam>_alpha_<alpha hyperparam>_<dd-mmm-yyyy>_<input file name>.mat`. We also have an example output of using `data.mat` as input in the repo named `beta_80_alpha_95%_31-Dec-2020_data.mat`.

Once you have the output file, you can load it into the Matlab and use some of our plotting function to visualize the tree. If there are less or equal to 2 TFs being used in the tree, you can use function `PlotTreeScatter`. In the input parameter, you will be asked if you want to plot each leaf as a single mean Gaussian or a linear regression. These are calculated using maximum likelihood after the Maximum A Posteriori tree is found by MCMC. Also, you can specify the number of bootstrap resampling to build an empirical confidence interval.

For example, if you just want each leaf to have a single mean:

```matlab
>> PlotTreeScatter(T_max, TF_level_log, gene_level_log, training_idx, TF_name, gene_name, Error, 0.95, 80, false, 100)
```
![Tree scatter plot with Gaussian leaf](/img/img2.png)

Or you want each leaf to have a linear regression:

```matlab
>> PlotTreeScatter(T_max, TF_level_log, gene_level_log, training_idx, TF_name, gene_name, Error, 0.95, 80, true, 100)
```

![Tree scatter plot with linear regression leaf](/img/img3.png)

Obviously when the tree uses more than 2 TFs, we cannot visualize it as in a 3D space. In this case you may use `PlotTreeStructure` function to visualize the tree:

```matlab
>> PlotTreeStructure(T_max, TF_level_log, gene_level_log, training_idx, TF_name, gene_name, median_factor, false)
```

![Tree violin plot](/img/img4.png)

Here each node has a violin plot of the distribution of the target gene. The black outline are the distribution of all data samples, and the colored part are the distribution of data samples that passed through the node. Each node has a title indicating the node ID, decision feature (TF name) and decision value. So in this plot, the mean of gene is 8.43 when TFE is higher than 5.82 and TFD lower than 4.61.

## Final Comments

There are a couple of unpleasant things about this model. Here we name two:

* Mixing is very slow. This is alluded to in the original paper, but in our case, only a few dozens of proposed trees are accepted in an MCMC run of a few thousand iterations. The acceptance rate is so low that we hardly think the MCMC is sampling the mode region of the posterior. Instead, it becomes more of an algorithm trying to **search for MAP** from the posterior.
* Although reducing a tree into a probability measure is nice for computation, we cannot really get a distribution of trees from the posterior. The posterior distribution, if we can find, represents trees of all kinds, with different number of nodes and edges. It would be easy to find the average(expectation) from the posterior distribution, but hard to find the "average" tree from the posterior. This is particularly unpleasant given we are working in Bayesian framework, although this makes the previous point less unpleasant. 
