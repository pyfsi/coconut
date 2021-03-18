# Models

This documentation describes the different types of available models. The purpose of a model is always to approximate a (inverse) Jacobian of a system, based on secant information from input and output pairs.
In order to approximate the Jacobian $\mathcal{A}'$ of a general function $a=\mathcal{A}(b)$, the model needs to be supplied with matching input and output pairs, ($b^i$, $a^i=\mathcal{A}(b^i)$).
Once at least two pairs have been supplied, the model is able to approximately predict the product of the Jacobian with an arbitrary vector $\Delta b$.
In other words, when the a vector $\Delta b$ is given it outputs $\Delta a=\widehat{\mathcal{A}}'\Delta b$, where the hat symbol is used to denote that an approximation of the Jacobain is used.

Which Jacobian is approximated in practice will depend on the use of the model in the coupled solver.
For example, using the coupled solver `CoupledSolverIQNI` with the model `ModelLS` corresponds to the IQN-ILS method developed by Degroote et al. [[1](#1)]. In that case, the approximated Jacobian is $\mathcal{R}'^{-1}$.
If the coupled solver `CoupledSolverIBQN` is combined with two instances the model `ModelLS`, the resulting algorithm corresponds to the IBQN-LS method developed by Vierendeels et al. [[2](#2)]. Then, the two models each approximate one Jacobian: $\mathcal{F}'$ and $\mathcal{S}'$.
Refer to the [coupled solvers documentation](../coupled_solvers.md) for more information on these notations.

In the following, the example from IQN-ILS will be used: the inverse Jacobian of $\mathcal{R}'$ with respect to $\tilde{x}$ is approximated which has an input vector $r$ and an output vector $\tilde{x}$.
For brevity, the approximation will denoted by $N^k$, where the superscript $k$ referes to the iteration.

## Common methods

There are three model-specific methods which are implemented by all models.

The first of which is the `predict(dr)` method, which returns an estimation of $\Delta \tilde{x}=N^k\Delta r$ from an input $\Delta r$, based on stored input and output pairs.
Second, in order to improve the estimation, input-output pairs can be added to the model using the method `add(r, xt)`.
Third, the method `filterq(dr)` returns the part of vector $\Delta r$ which is orthogonal to the columnspace of the matrix containing the differences between consecutively stored inputs.
In other words, it returns the part of the input vector $\Delta r$ for which the deficient approximation of the Jacobian holds no information. Note that it is equal to the part of the vector $\Delta r$ which is inside the nullspace of the approximation of the Jacobian represented by the model.

## Least-squares

The `type` for this model is `coupled_solvers.models.ls`.
The abbreviation LS stands for _least-squares_.
From the received input vectors $r^i$ and output vectors $\tilde{x}^i$, differences are constructed
$$
\delta r^{i-1}=r^i-r^{i-1},
\;
\delta \tilde{x}^{i-1}=\tilde{x}^i-\tilde{x}^{i-1}.
$$
These are stored as column of the matrices $V$ and $W$.
This model requires the approximation of the Jacobian, denoted by $N^k$, to fulfill the secant equations:
$$
W=N^k V.
$$
In addition to this a minimal norm requirement is imposed (hence the name).

The differences can be limited to the current time step only. However, using the reuse parameter `q`, the differences from the `q` previous time steps are included as well.
It is important to note that differences between the input and outputs from different time steps are not considered.
Reuse may greatly improve the convergence speed and stability. The optimal value of `q` is problem dependent. Typically, however, an optimal value is around 10.

This model is matrix-free, due to an implementation using QR-decomposition. With matrix-free is meant that no large dense matrices are constructed, not that no matrices are used at all.
The $R$ matrix from the QR-decomposition has to be invertible. Therefore, (almost) linearly dependent columns in the matrix containing the input information from the current and previous time steps should be removed. This is called filtering. The larger `q`, the more important filtering becomes.
If the diagonal element in $R$ is smaller than an absolute tolerance level `min_significant`, the corresponding column is removed.
The implementation is as such that the most recent information is kept.

For more information refer to [[3](#3)].

As mentioned before, the combination of this model wiht the coupled solver `CoupledSolverIQNI` corresponds to the IQN-ILS method developed by Degroote et al. [[1](#1)], while using twice this model with the coupled solver `CoupledSolverIBQN` corresponds to the IBQN-LS method developed by Vierendeels et al. [[2](#2)].

### Settings

The following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
<nobr>`min_significant`</nobr>|double|Absolute tolerance for filtering. To disable filtering set to `0`.
`q`|int|Number of previous time steps that are reused. In a steady simulation there are no previous time steps, so then it should be `0`.

## Multi-vector

The `type` for this model is `coupled_solvers.models.mv`. The abbreviation MV stands for _multi-vector_.

In this model, differences are constructed similar to the least-squares model. However, it requires the approximation $N^k$ to fulfill the secant equations of the current time step only. Moreover, it is required that the approximation is as close as possible to the previous time step. In this model large dense matrices are constructed and is hence discouraged for cases with a large number of degrees of freedom on the interface.

Filtering can also be applied here, then (almost) linear columns in the matrix containing the input information from the current time step are removed. However, filtering is much less critical compared to the `LS` model as it concerns only the information from the current time step. If no filtering is wanted, the tolerance level should be set to zero.

For more information refer to [[3](#3)].

The combination of this model with the coupled solver `CoupledSolverIQNI` corresponds to the IQN-MVJ from Lindner et al. [[4](#4)], while using twice this model with the coupled solver `CoupledSolverIBQN` corresponds to the MVQN method developed by Bogaers et al. [[5](#5)].

### Settings

The following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
<nobr>`min_significant`</nobr>|double|(optional) Default: `0` (disabled). Absolute tolerance for filtering.

## Multi-vector matrix-free

The `type` for this model is `coupled_solvers.models.mv_mf`. The abbreviation MV-MF stands for _multi-vector matrix-free_.

By combining the QR-approach from the least-squares model with the time step wise storage of secant information, a matrix-free implementation of the multi-vector approach is obtained quite naturally. This implementation was also thought of by Spenke et al. [[6](#6)].

In this approach, the contribution of each time step is grouped into a separate term, where the most recent time step has priority over the later ones. Therefore, a parameter `q` is used to denote how many time steps are reused. Setting this parameter very large, this model will act the same as `MV`, which reuses all time steps. Generally, the performance will increase when this parameter is chosen larger, but so will be the computational cost. Nonetheless, this cost is much smaller compared to the non-matrix-free multi-vector approach.

Filtering can be applied similar to the above described models, but this is typically not necessary.

### Settings

The following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
<nobr>`min_significant`</nobr>|double|(optional) Default: `0` (disabled). Absolute tolerance for filtering.
`q`|int|Number of previous time steps that are reused. In a steady simulation there are no previous time steps, so then it should be 0.

## References 
<a id="1">[1]</a> 
[Degroote J., Bathe K.-J., Vierendeels J., "Performance of a new partitioned procedure versus a monolithic procedure in fluid-structure interaction", Computers & Structures, vol. 87, no. 11–12, pp. 793-801, 2009.](http://hdl.handle.net/1854/LU-533365)

<a id="2">[2]</a> 
[Vierendeels J., Lanoye L., Degroote J., Verdonck P., "Implicit coupling of partitioned fluid-structure interaction problems with reduced order models", Computers & Structures, vol. 85, no. 11–14, pp. 970–976, 2007.](http://hdl.handle.net/1854/LU-409369)

<a id="3">[3]</a> 
[Delaissé N., Demeester T., Fauconnier D. and Degroote J., "Comparison of different quasi-Newton techniques for coupling of black box solvers", in ECCOMAS 2020, Proceedings, Paris, France, 2021.](http://hdl.handle.net/1854/LU-8685199)

<a id="4">[4]</a> 
[Lindner F., Mehl M., Scheufele K., Uekermann B., "A comparison of various quasi-Newton schemes for partitioned fluid-structure interaction", in: B. Schrefler, E. Oñate, M. Papadrakakis (Eds.), 6th International Conference on Computational 975 Methods for Coupled Problems in Science and Engineering, pp. 477–488, 2015.](https://www.researchgate.net/publication/277077208_A_Comparison_of_various_Quasi-Newton_Schemes_for_Partitioned_Fluid-Structure_Interaction)

<a id="5">[5]</a> 
[Bogaers A., Kok S., Reddy B., Franz T., "Quasi-Newton methods for implicit black-box FSI coupling", ComputerMethods in AppliedMechanics and Engineering, vol. 279, pp. 113–132, 2014.](https://doi.org/10.1016/j.cma.2014.06.033)

<a id="6">[6]</a> 
[Spenke T., Hosters N., Behr M., "A multi-vector interface quasi-newton method with linear complexity for partitioned fluid–structure interaction", Computer Methods in Applied Mechanics and Engineering, vol. 361, pp. 112810, 2020.](https://doi.org/10.1016/j.cma.2019.112810)