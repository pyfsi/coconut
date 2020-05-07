# Models

This documentation describes the models which are available.
These models always approximate a (inverse) Jacobian, denoted here by $J$.
Receiving an input, called $\Delta x$, they return an estimation of $J\Delta x=\Delta r$.
In order to improve the estimation, information is added, when it becomes available.
This information consists of pairs of calculated differences $\Delta x$ and $\Delta r$.
The different models use this information in a different way.

## LS

The abbreviation LS stands for _least-squares_.
This model requires the approximation of $J$, to fulfill the secant equations of the current and `q` previous time steps.
Moreover, it is required that the approximation is least-squares.
This model is matrix-free, due to an implementation using QR-decomposition.
With matrix-free is meant that no large dense matrices are constructed, not that no matrices are used at all.

The optimal value of `q` is problem dependent.
Typically however, an optimal value is around 10.
The $R$ matrix from the QR-decompostion has to be invertible.
Therefore, (almost) linear columns in the matrix containing the input information from the current and previous time steps should be removed.
This is called filtering. The larger `q`, the more important filtering becomes.
If the diagonal element in $R$ is smaller than an absolute tolerance level, the corresponding column is removed.
The most recent information is kept.

For more information refer to .

### Settings

The following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`min_significant`|double|Absolute tolerance for filtering.
`q`|int|Number of previous time steps that are reused. In a steady simulation there are no previous time steps, so then it should be 0.

## MV

The abbreviation MV stands for _multi-vector_.
not matrix free, not for large number of degrees of freedom on the interface.

This model requires the approximation of $J$, to fulfill the secant equations of the current time step only.
Moreover, it is required that the approximation is as close as possible to the previous time step.
In this model large dense matrices are constructed and is hence discouraged for cases with a large number of degrees of freedom on the interface.

Here filtering can also be applied.
Then, (almost) linear columns in the matrix containing the input information from the current time step are removed.
If the diagonal element in $R$ is smaller than an absolute tolerance level, the corresponding column is removed.
The most recent information is kept.
However, filtering is much less critical compared to the `LS` model as it concerns only the information from th current time step.
If no filtering is wanted, the tolerance level should be set to zero.

For more infromation refer to .

### Settings

The following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`min_significant`|double|Absolute tolerance for filtering.

## MV-MF

The abbreviation MV-MF stands for _multi-vector matrix-free_.
This method implements the multi-vector method, but in a matrix-free way.
Therefore, a parameter `q` is used to denote how many time steps are re-used.
Setting this parameter very large, this model will act the same as `MV`, which reuses all time steps.

### Settings

The following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`min_significant`|double|Absolute tolerance for filtering.
`q`|int|Number of previous time steps that are reused. In a steady simulation there are no previous time steps, so then it should be 0.
