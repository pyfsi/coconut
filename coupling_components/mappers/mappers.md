# Mappers



The coupling interfaces in different solvers are often discretized (meshed) differently, resulting in non-conforming `ModelParts`. For this reason, mapping is usually required, to transfer data from the output `Interface` of one solver to the input `Interface` of another. 

The terms *from* and *to* will be used to denote respectively the side that provides the original data and the side that receives the mapped data.  


## General concepts


### Hierarchy of mapping-related objects

CoCoNuT interacts with the mappers through an instance of the `SolverWrapperMapped` class. This special solver wrapper appears and behaves exactly the same as real solver wrappers.
It contains 3 `Components`: a mapper for the input, a real solver wrapper and a mapper for the output.

The two mappers in the `SolverWrapperMapped` are of a special type: they work on the level of the [`Interfaces`](../../data_structure/data_structure.md#interface), while the actual mapping will be done on a lower level: the [`ModelPart`](../../data_structure/data_structure.md#modelpart) level. The `Interface` mapper effectively acts as a *wrapper* around the actual `ModelPart` mappers. Currently only one mapper is available on the `Interface` level, aptly called `MapperInterface`.

At the lowest level, mappers interpolate variable data (stored in `Interfaces`) between two `ModelParts`, based on their coordinates. Interpolation is always done from the *from* `ModelPart` to the *to*-`ModelPart`.
Multiple `ModelPart` mappers can be chained together in a `MapperCombined` object, to add additional functionality, for example swapping the coordinate axes with the [permutation mapper](mappers.md#mapperpermutation).


### Interpolators and transformers

The `ModelPart`-level mappers have two main methods: `initialize` and `__call__`. 
The `initialize` method performs one-time expensive operations, such as a nearest-neighbour search and the calculation of interpolation coefficients. 
The `__call__` method is used for the actual mapping. It takes two tuples as arguments (*from* and *to* respectively). Each tuple contains an `Interface` object, the name of the affected `ModelPart` and the name of the variable that must be mapped. This method returns nothing: the mapping is done in-place in the *to*-`Interface`.

There are two types of `ModelPart` mappers: _interpolators_ and _transformers_. They can be distinguished by their superclass. 
All interpolators inherit from the superclass `MappedInterpolator`. These mappers do the actual interpolation. Currently [a nearest-neighbour mapper](mappers.md#mappernearest), [a linear mapper](mappers.md#mapperlinear) and a [radial basis mapper](mappers.md#mapperradialbasis) are available.
All transformers inherit from the superclass `MappedTransformer`. They provide additional functionality that can be useful during mapping. They do not map values from one point cloud to another like the interpolators, but perform *one-sided* transformations on the coordinates and/or the data. One example is the permutation transformer: it swaps the coordinate axes of the `ModelPart` and accordingly the components of vector variables. 
A transformer can never be used by itself, it must always be combined with an interpolator. The reason is that interpolators use information that comes from two sides, which is exactly what the higher-level `SolverWrapperMapped` and `MapperInterface` objects are supplying. To chain together multiple `ModelPart` mappers, the `MapperCombined` is used: it contains always one interpolator and zero or more transformers, on either side of the interpolator.



## Special mapper classes


### `MapperInterface`

Class that maps on the level of `Interface` objects. 
It takes two `Interfaces`, and maps the `ModelParts` to each other in the order in which the `Interfaces` are defined in the JSON file.
The same type of `ModelPart` mapper is used for each pair of `ModelParts`: this can either be an interpolator or a combined mapper. To use a different `ModelPart` mapper for the different `ModelParts` in the `Interface`, or even for different variables, a new `Interface` mapper would have to be written. 

parameter|type|description
------:|:----:|-----------
`type`|str|`ModelPart` mapper to be used.
<nobr>`settings`</nobr>|dict|All the settings for the `ModelPart` mapper specified in `type`.


### `MapperInterpolator`

Superclass for all [interpolators](mappers.md#interpolators). 

parameter|type|description
------:|:----:|-----------
`directions`|list|List of coordinate directions, maximum three entries, may contain `"x"`, `"y"`, `"z"`.
`scaling`|list|Optional. List of scaling factors, must be same length as `directions`. Coordinates are scaled with these factors, this may improve interpolation e.g. when cells have a high aspect ratio with respect to one of the axes. 
<nobr>`balanced_tree`</nobr>|bool|Optional, default `false`. If set to `true` a balanced `cKDTree` is created, which is more stable, but takes longer to build. Set to `true` in the rare case that the tree gives problems.

The `initialize` method must be defined in all child classes. It takes as arguments the *from*-`ModelPart` and the *to*-`ModelPart`. It does the following:

-   read and store the coordinates from the *from*- and *to*-`ModelParts`,
-   scale coordinates if necessary,
-   check if the bounding boxes of the *from*- and *to*-`ModelParts` overlap,
-   do an efficient nearest neighbour search using `scipy.spatial.cKDTree`,
-   check if the *from*-`ModelPart` does not contain duplicate coordinates.

The `__call__` method should not be overridden in the child classes. It interpolates data based on coefficients that are calculated in the `initialize` method of the child classes. Both scalar and vector variables can be mapped.



###  `MapperTransformer`

Superclass for all [transformers](mappers.md#transformers). A transformer cannot be used as a standalone mapper, but will always be part of a combined mapper. 

The `initialization` of transformers is very different from that of interpolators. A transformer is initialized from one side (either the *from* or the *to* side), i.e. based on a single `ModelPart`. From this input `ModelPart`, the `initialize` method creates and returns another `ModelPart` that is stored and used in the combined mapper. 



### `MapperCombined`

parameter|type|description
------:|:----:|-----------
<nobr>`mappers`</nobr>|list|An ordered list of all the `ModelPart` mappers to be used.

The `MapperCombined` is used to chain together multiple mappers. It always contains a single interpolator and zero or more transformers, which can be on either side of the interpolator. If transformers are present, _intermediate_ `ModelParts` are created during initialization. This is done by working _inwards_ towards the interpolator. This means that transformers upstream of the interpolator, are initialized based on the *from*-`ModelPart` (input), while downstream transformers are initialized based on the *to*-`ModelPart` (output). 
Some transformers can only be initialized in one direction, e.g. for `MapperAxisymmetric3DTo2D`, the 2D *to*-`ModelPart` must be supplied, therefore this transformer must be downstream of the interpolator. 

These concepts are clarified further using an excerpt from the JSON file of the *Fluent 3D - Abaqus 2D* tube example: 

```JSON
{
"type": "mappers.combined",
"settings": {
    "mappers": [
        {"type": "mappers.permutation",
        "settings": {"permutation": [1, 0, 2]}},
        {"type": "mappers.radial_basis",
        "settings": {"directions": ["x", "y", "z"]}},
        {"type": "mappers.axisymmetric_3d*to*2d",
        "settings": {"direction_axial": "y", "direction_radial": "x", "n_tangential": 8}}
        ]
    }
}
```
This combined mapper contains 3 `ModelPart` mappers. In the `initialize` method, the following happens (in this order):

-   A first intermediate `ModelPart` is created by calling the `initialize` method of the permutation transformer. This is done by swapping the x and y coordinates of the *from*-`ModelPart`. 
-   A second intermediate `ModelPart` is created by calling the `initialize` method of the axisymmetric transformer. The coordinates of the 2D *to*-`ModelPart` are used to create the 3D intermediate `ModelPart` by adding points in the circumferential direction.
-   The mapping coefficients of the radial basis interpolator are calculated by calling its `initialize` method. The first intermediate `ModelPart` is used as *from*-`ModelPart`, the second intermediate `ModelPart` is used as *to*-`ModelPart`.

When the `__call__` method  of the combined mapper is used, the following happens (in this order):

-   The *from*-data (stored in the *from*-`Interface`) is mapped from the *from*-`ModelPart` to the first intermediate `ModelPart` using the `__call__` method of the permutation mapper: scalar variables are unchanged, vector variables are permuted. 
-   The resulting data is now interpolated to the second intermediate `ModelPart`, using the `__call__` method of the radial basis mapper. 
-   Finally that data is mapped to the *to*-`ModelPart` using the `__call__` method of the axisymmetric transformer, reducing it from 3D to 2D. That data is written to the *to*-`Interface`. 





## Transformers


### `MapperPermutation`

Permutates the coordinates and the vector variables according to the given `permutation`. 
This transformer can be initialized in both directions. 

parameter|type|description
------:|:----:|-----------
<nobr>`permutation`</nobr>|list|A permutation of the list [0, 1, 2].


### `MapperAxisymmetric2DTo3D`

Transforms a 2D axisymmetric geometry to a 3D geometry. This transformer can only be initialized in the _forward_ direction, i.e. based on the 2D axisymmetric `ModelPart`. Therefore, it should be _upstream_ of the interpolator in the combined mapper.

The 3D `ModelPart` is returned by the initialization.`angle` defines the 3D geometry and is based on the geometry of the 3D solver, e.g. the geometry of an axisymmtric simulation with an OpenFoam solver is defined with an angle of 5°.  `n_tangential` specifies the number of points in the tangential (circumferential) direction, so that for each point in the 2D `ModelPart`, `n_tangential` points are created in the 3D `ModelPart`. It is important to make the value of `n_tangential` large enough, ideally close to the number of points in the tangential direction in the 3D solver. The code knows which directions are axial, radial and tangential thanks to the input parameters `direction_axial` and `direction_radial`. 

It is not possible to change the axial direction between 2D and 3D: a separate `MapperPermutation` should be added to the combined mapper for that purpose. 

Scalar data is simply copied from the 2D point to all corresponding 3D points. For vector data, the axial component is simply copied, the radial component is rotated. The tangential component (e.g. swirl) is not taken into account currently.

Points that lie on the symmetry axis can not be handled by the current transformer.


parameter|type|description
------:|:----:|-----------
`direction_axial`|str|Must be `"x"`, `"y"` or `"z"`, specifies the symmetry axis.
<nobr>`direction_radial`</nobr>|str|Must be `"x"`, `"y"` or `"z"`, specifies the second (radial) axis in 2D.
`n_tangential`|int|Degrees of freedom in tangential (circumferential) direction of 3D `ModelPart` that is created during initialization. The minimum setting of n_tangential points depends of the definition of the angle.
`angle`|int|(optional) Default: `360`. Angle of the (partial) 3D cylinder constructed from the 2D geometry, centred around the radial direction.  


### `MapperAxisymmetric3DTo2D`

Transforms a 3D geometry to a 2D axisymmetric geometry. This transformer can only be initialized in the _backward_ direction, i.e. based on the 2D axisymmetric `ModelPart`. Therefore, it should be _downstream_ of the interpolator in the combined mapper.

For scalar data, the circumferential average is taken for each 2D point. For vector data too, taking into account the correct radial direction in each 3D point.  Again, swirl cannot be taken into account: if a tangential component is present in 3D, it is not transferred to 2D. 

More information and JSON settings can be found under [`MapperAxisymmetric2DTo3D`](mappers.md#mapperaxisymmetric2dto3d).


### `MapperDepth2DTo3D`

Transforms a 2D geometry to a 3D geometry, by extruding in a _depth_ direction. This transformer can only be initialized in the _forward_ direction, i.e. based on the 2D axisymmetric `ModelPart`. Therefore, it should be _upstream_ of the interpolator in the combined mapper.

The 3D `ModelPart` is returned by the initialization. The depth direction in which the 2D `ModelPart` is extended, is specified by `direction_depth`. This direction has to be along one of the principal axes: x, y and z. The number of nodes in the depth direction and their location are given by the list `coordinates_depth`.

Scalar data is simply copied from the 2D point to all corresponding 3D points. For vector data, all components other than the depth component are simply copied. The depth component is set to zero.


parameter|type|description
------:|:----:|-----------
<nobr>`coordinates_depth`</nobr>|list|Contains the depth coordinates to which the nodes of the orignal 2D plane are copied.
`direction_depth`|str|Must be `"x"`, `"y"` or `"z"`, specifies the symmetry axis.


### `MapperDepth3DTo2D`

Transforms a 3D geometry to a 2D axisymmetric geometry, by collapsing in a _depth_ direction. This transformer can only be initialized in the _backward_ direction, i.e. based on the 2D axisymmetric `ModelPart`. Therefore, it should be _downstream_ of the interpolator in the combined mapper.

For scalar data, the average is taken for each 2D point. For vector data too, again removing the depth component.

More information and JSON settings can be found under [`MapperAxisymmetric2DTo3D`](mappers.md#mapperaxisymmetric2dto3d).


## Interpolators



### `MapperNearest`

Does not require additional settings compared to the `MapperInterpolator`. Does simple nearest-neighbour mapping.


### `MapperLinear`

Additional settings:

parameter|type|description
------:|:----:|-----------
<nobr>`parallel`</nobr>|bool|Optional, default `false`. If `true` the package `multiprocessing` is used to parallellize the loop that the calculates the interpolation coefficients. This is only useful for `ModelParts` with a very high number of degrees of freedom. 

The kind of linear mapping depends on the number of coordinate directions, as given in the `directions` setting.

**1D:** If the *to*-point lies between the 2 nearest *from*-points, linear interpolation is done. Else, nearest neighbour interpolation is done.

**2D:** The *to*-point is first projected on the line through the 2 nearest *from*-points. If the projected point lies between the *from*-points, linear interpolation is done. Else, nearest neighbour interpolation is done.

**3D:** The *to*-point is first projected on the plane through the 3 nearest *from*-points. If the triangle consisting of those 3 points is _deprecated_ (colinear points), the 2D-methodology is followed. Else, if the projected point lies inside the triangle, barycentric interpolation is done. If it lies outside the triangle, the 2D-methodology is followed.


### `MapperRadialBasis`

Additional settings:

parameter|type|description
------:|:----:|-----------
<nobr>`parallel`</nobr>|bool|Optional, default `false`. If `true` the package `multiprocessing` is used to parallellize the loop that the calculates the interpolation coefficients. This is only useful for `ModelParts` with a very high number of degrees of freedom.
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
 
Assume that $n$ nearest *from*-points will be used in the interpolation.
An unknown function $f(\boldsymbol{x})$ can then be approximated as the weighted sum of $n$ shifted radial basis functions:

$$
f(\boldsymbol{x}) = \sum_j \alpha_j \phi(||\boldsymbol{x} - \boldsymbol{x}_j||)
$$

To determine the coefficients $\alpha_j$, we require that the exact function value is returned at the $n$ *from*-points.
This gives us $n$ equations

$$
f(\boldsymbol{x}_i) = f_i = \sum_j \alpha_j \phi(||\boldsymbol{x}_i - \boldsymbol{x}_j||)
$$

which can be written in matrix form as

$$
\boldsymbol{f} = \boldsymbol{\Phi} \boldsymbol{\alpha}
$$

with $\boldsymbol{f}, \boldsymbol{\alpha} \in \mathbb{R}^{n \times 1}$, and $\boldsymbol{\Phi} \in \mathbb{R}^{n \times n}$. This system can be solved for the weights-vector $\boldsymbol{\alpha}$.

However, in our case, the *from*-point values vector $\boldsymbol{f}$ is not known in advance: it contains the values of the `Variable` that will be interpolated. 

Therefore, the approximation to calculate the interpolatoin in the *to*-point is rewritten as follows:

$$
f(\boldsymbol{x}_{to}) = \sum_j \alpha_j \phi(||\boldsymbol{x}_{to} - \boldsymbol{x}_j||) = \boldsymbol{\Phi}_{to}^T \boldsymbol{\alpha} = \boldsymbol{\Phi}_{to}^T \boldsymbol{\Phi}^{-1} \boldsymbol{f} = \boldsymbol{c}^T \boldsymbol{f}
$$

The coefficients vector $\boldsymbol{c}$ can now be calculated based only on the coordinates by solving the system

$$
\boldsymbol{\Phi} \ \boldsymbol{c} = \boldsymbol{\Phi}_{to}.
$$

As every to-point has different nearest neighbours in the *from*-points, the coefficient vector $\boldsymbol{c}$ must be calculated for each *to*-point independently. The matrix $\boldsymbol{\Phi}$ and vector $\boldsymbol{\Phi}_{to}$ must also be calculated for every *to*-point independently.

For every *to*-point, the reference distance $d_{ref}$ is determined as the product of the `shape_parameter` and the distance between the *to*-point and the furthest *from*-point.

In order to ensure that the basis function of each of the nearest *from*-points covers every *from*-point, the `shape_parameter` should be larger than two.
This value may however lead to an interpolation function which consists of sharp peaks or wiggles, with the correct value near the *from*-points, but a deviating value away from them.

In the extreme case of $d_{ref}$ approaching zero, the so-called "bed-of-nails interpolant" is obtained, which is close to zero everywhere, except near the *from*-points where it sharply peaks.
In this case the interpolation matrix approaches the identity matrix.

Choosing a higher value improves the interpolation as the basis functions become wider, but the interpolation matrix becomes less stable, i.e. the condition number increases.
The default value is 200.
In practice, the `shape_parameter` is chosen so that the interpolation matrix is "on the edge of ill-conditioning" (eg. with a condition number of roughly $10^{13}$ for double-precision floating point).
A warning is printed when the condition number of an interpolation matrix becomes higher than $10^{13}$.
