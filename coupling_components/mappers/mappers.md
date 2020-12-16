# Mappers


## General concepts


### Hierarchy of mapping-related objects

CoCoNuT interacts with the mappers through the `SolverWrapperMapped` object. This solver-wrapper appears and behaves exactly the same as real solver-wrappers.
It contains 3 `Components`: a mapper for the input, a real solver-wrapper and a mapper for the output. The mappers are initialized through the `set_interface_input` and `set_interface_output` methods respectively, by providing them with the `Interfaces` that will be respectively the input and output of the `SolverWrapperMapped` object.

The two mappers in the `SolverWrapperMapped` are also of a special type: they work on the level of `Interfaces`. They are some sort of _mapper-wrapper_ around the actual mappers which work on `ModelPart` level.
Currently only one such mapper is available, aptly called `MapperInterface`.

At the lowest level, mappers interpolate historical variables between two `ModelParts`, based on the coordinates of the `Nodes`. Interpolation is always done from the _from_-`ModelPart` to the _to_-`ModelPart`.
These mappers can be chained together in a `MapperCombined` mapper, creating in fact another layer of mapping. So many layers! Like an onion!


### Interpolators and transformers

The `ModelPart`-level mappers have two main methods: `Initialize` and `__call__`. 

The `Initialize` method performs one-time expensive operations, namely nearest-neighbour search and calculation of the interpolation coefficients. The initialization is done based on the original coordinates `X0`, `Y0` and `Z0`.

The `__call__` method is used for the actual mapping. It takes two tuples as arguments (_from_ and _to_ respectively), each tuple containing the `ModelPart` and the `Variable` to be used in the interpolation. This method returns nothing: the interpolation is done in-place in the `ModelParts`.

There are two types of `ModelPart`-level mappers: _interpolators_ and _transformers_. They can be distinguished by their boolean `interpolator` attribute (see `__init__`). 

For interpolators, `Initialize` gets two `ModelPart` objects (_from_ and _to_), and returns nothing. These mappers do the real interpolation. Currently `MapperNearest`, `MapperLinear` and `MapperRadialBasis` are available.

For transformers, `Initialize` gets only one `ModelPart` (from either the _from_ or _to_ side, depending on the transformation), and returns the other `ModelPart`. Currently `MapperPermutation`, `MapperAxisymmetric2DTo3D` and `MapperAxisymmetric3DTo2D` are available.

A transformer can never be used by itself, it must always be combined with an interpolator: the reason is that interpolators use information coming from two sides, which is exactly what the `SolverWrapperMapped` and `MapperInterface` objects want. To chain together multiple mappers, the `MapperCombined` is used: it contains always 1 interpolator and 0 or more transformers, on either side of the interpolator.


## Overview of special mappers


### MapperInterface

Special mapper-class that maps on the level of `Interface` objects. 
It takes two `Interfaces`, and maps the `ModelParts` to each other in order of appearance, all using the same `ModelPart` mapper.

To use a different `ModelPart` mapper for the different `ModelParts` in the `Interface`, or even for the different historical variables, a new `Interface` mapper must be written. 

JSON setting|type|description
------:|:----:|-----------
`type`|str|`ModelPart` mapper to be used.
<nobr>`settings`</nobr>|dict|All the settings for the `ModelPart` mapper specified in `type`.



### MapperCombined

The `MapperCombined` is used to chain together multiple mappers. It contains always 1 interpolator and 0 or more transformers, on either side of the interpolator. If transformers are present, _intermediate_ `ModelParts` are created during initialization. This is done by working _inwards_ towards the interpolator. This means that transformers upstream of the interpolator, are initialized based on the _from_ `ModelPart` (input), while downstream transformers are initialized based on the _to_ `ModelPart` (output). 
Some transformers can only be initialized in one direction, e.g. for `MapperAxisymmetric3DTo2D`, the 2D _to_-`ModelPart` must be supplied, therefore it must be downstream of the interpolator. 

JSON setting|type|description
------:|:----:|-----------
<nobr>`mappers`</nobr>|list|An ordered list of all the `ModelPart` mappers to be used.



## Overview of transformers

### MapperPermutation

Permutates the coordinates and the vector variables according to the given `permutation`. 
This transformer can be initialized in both directions. 

JSON setting|type|description
------:|:----:|-----------
<nobr>`permutation`</nobr>|list|A permutation of the list [0, 1, 2].

### MapperAxisymmetric2DTo3D

Transforms from a 2D axisymmetric geometry to a 3D geometry. This transformer can only be initialized in the _forward_ direction, i.e. based on the 2D axisymmetric `ModelPart`. Therefore, it should be _upstream_ of the interpolator in the `MapperCombined`.

The 3D `ModelPart` is returned by the initialization. `n_tangential` specifies the number of `Node` objects in the tangential (circumferential) direction. For each `Node` in the 2D `ModelPart`, `n_tangential` `Nodes` are created in the 3D `ModelPart`. It is important to take the value of `n_tangential` large enough. By preference as close as possible to the number of nodes in the tangential direction of the 3D `ModelPart`. The code knows which directions are axial, radial and tangential thanks to the input parameters `direction_axial` and `direction_radial`. 

It is not possible to change the axial direction between 2D and 3D: a separate `MapperPermutation` should be added for that purpose. 

Scalar variables are simply mapped from the 2D `Node` to all corresponding 3D `Nodes`. For vector variables, the axial component is simply mapped, the radial component is rotated. The tangential component (e.g. swirl) cannot be taken into account.


JSON setting|type|description
------:|:----:|-----------
`direction_axial`|string|Must be `"X"`, `"Y"` or `"Z"`, specifies the symmetry axis.
<nobr>`direction_radial`</nobr>|string|Must be `"X"`, `"Y"` or `"Z"`, specifies the second (radial) axis in 2D.
`n_tangential`|int|Degrees of freedom in tangential (circumferential) direction of 3D `ModelPart` that is created during initialization. Must be $\geq 6$.

### MapperAxisymmetric3DTo2D

Transforms from a 3D geometry to a 2D axisymmetric geometry. This transformer can only be initialized in the _backward_ direction, i.e. based on the 2D axisymmetric `ModelPart`. Therefore, it should be _downstream_ of the `MapperInterpolator` in the `MapperCombined`.

For scalar variables, the circumferential average is taken for each 2D `Node`. For vector variables too, taking into account the correct radial direction in each 3D `Node`.  Again, swirl cannot be taken into account: if a tangential component is present in 3D, it is not mapped to 2D. 

For more information and JSON settings, see `MapperAxisymmetric2DTo3D` which is very similar.



## Overview of interpolators

### MapperInterpolator

Base-class for all interpolators (currently `MapperNearest`, `MapperLinear` and `MapperRadialBasis`). 

JSON setting|type|description
------:|:----:|-----------
`directions`|list|List of coordinate directions, maximum three entries, may contain `"X"`, `"Y"`, `"Z"`.
`scaling`|list|Optional. List of scaling factors, must be same length as `directions`. Coordinates are scaled with these factors, this may improve interpolation e.g. when cells have a high aspect ratio with respect to one of the axes. 
<nobr>`balanced_tree`</nobr>|bool|Optional, default `false`. If set to `true` a balanced `cKDTree` is created, which is more stable, but takes longer to build. Set to `true` in the rare case that the tree gives problems.

The `Initialize` method should be called in all child-classes. It does the following:

-   read and store the coordinates from the _from_ and _to_ `ModelParts`
-   scale coordinates if necessary
-   check if the bounding boxes of the _from_ and _to_ `ModelParts` are more or less overlapping
-   do an efficient nearest neighbour search using `scipy.spatial.cKDTree`
-   check if the _from_ `ModelPart` does not contain duplicate `Nodes` (i.e. with same coordinates)

The `__call__` method should not be overridden in the child-classes. It maps historical variables based on neighbours and coefficients determined in `Initialize`. Historical variables of type `Double` and type `Array` can be mapped (the latter is just the application of the former for each vector component).


### MapperNearest

Child-class of `MapperInterpolator`, does not require additional settings. Does simple nearest-neighbour mapping.


### MapperLinear

Child-class of `MapperInterpolator`, additional settings:

JSON setting|type|description
------:|:----:|-----------
<nobr>`parallel`</nobr>|bool|Optional, default `false`. If `true` the package `multiprocessing` is used to parallellize the loop that the calculates the interpolation coefficients. This is only useful for `ModelParts` with a very high number of `Nodes`. 

The kind of linear mapping depends on the number of coordinate directions, as given in the `directions` setting.

**1D** - If the _to_-point lies between the 2 nearest _from_-points, linear interpolation is done. Else, nearest neighbour interpolation is done.

**2D** - The _to_-point is first projected on the line through the 2 nearest _from_-points. If the projected point lies between the _from_-points, linear interpolation is done. Else, nearest neighbour interpolation is done.

**3D** - The _to_-point is first projected on the plane through the 3 nearest _from_-points. If the triangle consinsting of those 3 points is _deprecated_ (colinear points), the 2D-methodology is followed. Else, if the projected point lies inside the triangle, barycentric interpolation is done. If it lies outside the triangle, the 2D-methodology is followed.


### MapperRadialBasis

Child-class of `MapperInterpolator`, additional settings:

JSON setting|type|description
------:|:----:|-----------
<nobr>`parallel`</nobr>|bool|Optional, default `false`. If `true` the package `multiprocessing` is used to parallellize the loop that the calculates the interpolation coefficients. This is only useful for `ModelParts` with a very high number of `Nodes`.
<nobr>`shape_parameter`</nobr>|int|Optional, default `200`. Should be chosen as large as possible without rendering the interpolation matrix ill-conditioned.

Radial basis function interpolation is relatively straightforward: implementation for 1D, 2D and 3D is exactly the same and can be written in a condensed way using `scipy.spatial.distance`. 

Normal radial basis interpolation is done as follows.
$\phi(r)$ is a radial basis function defined as  


$$
\phi(r) = 
\begin{cases}
    (1-\frac{r}{d_{ref}})^4 (1 + 4\frac{r}{d_{ref}}) \quad &\mathrm{for} \quad 0 \leq \frac{r}{d_{ref}} < 1 \\
    0 &\mathrm{for} \quad 1 \leq \frac{r}{d_{ref}}
\end{cases}
$$

with $r$ a positive distance. To control the width of the function, $r$ is scaled with a reference distance $d_{ref}$.
 
Assume that $n$ nearest _from_-points will be used in the interpolation.
An unknown function $f(\boldsymbol{x})$ can then be approximated as the weighted sum of $n$ shifted radial basis functions:

$$
f(\boldsymbol{x}) = \sum_j \alpha_j \phi(||\boldsymbol{x} - \boldsymbol{x}_j||)
$$

To determine the coefficients $\alpha_j$, we require that the exact function value is returned at the $n$ _from_-points.
This gives us $n$ equations

$$
f(\boldsymbol{x}_i) = f_i = \sum_j \alpha_j \phi(||\boldsymbol{x}_i - \boldsymbol{x}_j||)
$$

which can be written in matrix form as

$$
\boldsymbol{f} = \boldsymbol{\Phi} \boldsymbol{\alpha}
$$

with $\boldsymbol{f}, \boldsymbol{\alpha} \in \mathbb{R}^{n \times 1}$, and $\boldsymbol{\Phi} \in \mathbb{R}^{n \times n}$. This system can be solved for the weights-vector $\boldsymbol{\alpha}$.

However, in our case, the _from_-point values vector $\boldsymbol{f}$ is not known in advance: it contains the values of the `Variable` that will be interpolated. 

Therefore, the approximation to calculate the interpolatoin in the _to_-point is rewritten as follows:

$$
f(\boldsymbol{x}_{to}) = \sum_j \alpha_j \phi(||\boldsymbol{x}_{to} - \boldsymbol{x}_j||) = \boldsymbol{\Phi}_{to}^T \boldsymbol{\alpha} = \boldsymbol{\Phi}_{to}^T \boldsymbol{\Phi}^{-1} \boldsymbol{f} = \boldsymbol{c}^T \boldsymbol{f}
$$

The coefficients vector $\boldsymbol{c}$ can now be calculated based only on the coordinates by solving the system

$$
\boldsymbol{\Phi} \ \boldsymbol{c} = \boldsymbol{\Phi}_{to}.
$$

As every to-point has different nearest neighbours in the _from_-points, the coefficient vector $\boldsymbol{c}$ must be calculated for each _to_-point independently. The matrix $\boldsymbol{\Phi}$ and vector $\boldsymbol{\Phi}_{to}$ must also be calculated for every _to_-point independently.

For every _to_-point, the reference distance $d_{ref}$ is determined as the product of the `shape_parameter` and the distance between the _to_-point and the furthest _from_-point.

In order to ensure that the basis function of each of the nearest _from_-points covers every _from_-point, the `shape_parameter` should be larger than two.
This value may however lead to an interpolation function which consists of sharp peaks or wiggles, with the correct value near the _from_-points, but a deviating value away from them.

In the extreme case of $d_{ref}$ approaching zero, the so-called "bed-of-nails interpolant" is obtained, which is close to zero everywhere, except near the _from_-points where it sharply peaks.
In this case the interpolation matrix approaches the identity matrix.

Choosing a higher value improves the interpolation as the basis functions become wider, but the interpolation matrix becomes less stable, i.e. the condition number increases.
The default value is 200.
In practice, the `shape_parameter` is chosen so that the interpolation matrix is "on the edge of ill-conditioning" (eg. with a condition number of roughly $10^{13}$ for double-precision floating point).
A warning is printed when the condition number of an interpolation matrix becomes higher than $10^{13}$.
