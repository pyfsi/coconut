# Coupled solvers

This documentation describes the different types of coupled solvers which are available.
A coupled solver refers to a coupling algortihm used to couple two solvers, typically a flow and a structure solver.
Some of these coupled solvers make use of one or more [models](models/models.md).
An odd one out is `test_single_solver` which allows to test only one solver by coupling it to a dummy solver. 
All coupled solvers inherit from the class `GaussSeidel`.

In the parameter JSON file, the dictionary `coupled_solver` holds the `type` and the dictionary `settings`,
but also the dictionary `predictor`, the dictionary `convergence_criterion` and the list `solver_wrappers` containing 2 dictionaries: one for each `solver_wrapper`.
More information on these last three can be found in [predictors](../predictors/predictors.md),
[convergence criteria](../convergence_criteria/convergence_criteria.md) and the [solver_wrappers](../solver_wrappers/solver_wrappers.md) documentation, respectively.

In the following subsections, basics scheme will be shown. In those schematics, $\mathcal{F}$ is the first solver with input $x$ and output $\tilde{y}$ 
and $\mathcal{S}$ is the second solver with input $y$ and output $\tilde{x}$.
Typically, these solvers are a flow and structure solver, respectively.
The superscript $k=0\dots$ denotes the iteration, where $k+1$ is the current iteration.

## Gauss-Seidel

The `type` for this coupled solver is `coupled_solvers.gauss_seidel`.

### Algorithm

Gauss-Seidel or fixed-point iterations are the simplest way of coupling two solvers: the output of one solver is given to the other one without adjustment.

The following figure shows the basic methodology.
![gauss_seidel](images/iteration_gauss_seidel.png "Gauss-Seidel iterations")

Gauss-Seidel iterations are very simple, but unstable for cases with incompressible flow and high added mass.
A considerable convergence stabilization and acceleration is obtained by modifying the input to the one or both of the solvers using derivative information, as will be shown further.

### Settings

The following parameters need to be included in the `settings` dictionary.
Here they are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`delta_t`|double|Fixed time step size used in both solvers. For a steady simulation typically a value of 1 is taken.
`name`|string|(optional) Default: "results". Name of the case used to store a pickle file with results.
`save_results`|bool|(optional) Default: `false`. If `true` a [pickle](https://docs.python.org/3.9/library/pickle.html) file is stored containing some main [results](#save-results).
`time_step_start`|int|Time step number to (re)start a transient FSI calculation. If 0 is given, the simulation starts from scratch. Otherwise, the code looks for the relevant files to start from the corresponding time step. Not every solver wrapper implements restart, see the corresponding documentation for more information. For a steady simulation, the value should be 0.

`timestep_start` and `delta_t` are necessary parameters (also in a steady simulation), but can also defined in the solver wrapper directly (e.g. for standalone testing).
If they are defined both here and in the solver wrapper, then the former value is used and a warning is printed.

**These parameters should also be specified for the coupled solvers inheriting from `GaussSeidel`.**

## Relaxation

This coupled solver inherits from the class `GaussSeidel`.
The `type` for this coupled solver is `coupled_solvers.relaxation`.

### Algorithm

Gauss-Seidel iterations are very simple, but are unstable for cases with incompressible flow and high added mass.
A simple approach to mitigate this is applying relaxation, also called simple mixing.
In this coupling method the output of first solver is still given to the other one without adjustment, 
but the output of the second solver is relaxed as follows:
$$
x^{k+1}=(1-\omega)x^k+\omega\tilde{x}^k=x^k+\omega r^k
$$
with $x^k$ and $\tilde{x}^k$, respectively the input for the first solver and the output of the second solver in iteration $k$.
The difference between both is called the residual $r^k=\tilde{x}^k-x^k$.
The mixing or relaxation factor is $\omega$.

A symbolic schematic is given in the following figure.
![relaxation](images/iteration_relaxation.png "Relaxed Gauss-Seidel iterations")

This method is again quite simple, but be able to stabilize some cases which fail with Gauss-Seidel iterations.
A lower $\omega$ increases stability, but decreases convergence speed.
For more challenging problems, with incompresible flow and high added-mass, this approach will result in a very slow convergence, if it converges at all.

### Settings

Besides the parameters required in the class `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`omega`|double|Relaxation factor.

## Aitken

This coupled solver inherits from the class `GaussSeidel`.
The `type` for this coupled solver is `coupled_solvers.aitken`.

### Algorithm

Gauss-Seidel iterations are very simple, but are unstable for many cases.
An approach to mitigate this is applying relaxation.
However, the choice of a relaxation factor can be difficult: a high factor leads to a higher degree of stability, but a lower convergence speed.
In more challenging cases there will be no acceptable value for $\omega$.
In this coupling method a dynamic relaxation factor is used.
The output of the first solver is still given to the other one without adjustment, 
but the output of the second solver is relaxed as follows:
$$
x^{k+1}=(1-\omega^k)x^k+\omega^k\tilde{x}^k=x^k+\omega^k r^k
$$
with $x^k$ and $\tilde{x}^k$, respectively the input for the first solver and the output of the second solver in iteration $k$.
The difference between both is called the residual $r^k=\tilde{x}^k-x^k$.
The mixing or relaxation factor is $\omega^k$ is dynamic in the sense that its value changes between iterations.

A symbolic schematic is given in the following figure.
![aitken](images/iteration_aitken.png "Aitken iterations")

The value of $\omega^k$ is determined by applying the secant method for scalars directly to vectors and projecting it on $r^k-r^{k-1}$
as follows
$$
\omega^k=-\omega^{k-1}\frac{(r^{k-1})^T(r^k-r^{k-1})}{(r^k-r^{k-1})^T(r^k-r^{k-1})}.
$$

The first relaxation factor in a time step is equal to the last relaxation factor from the previous time step $\omega^n$, but limited to $\omega^{max}$.
$$
\omega^0=\textrm{sign}(\omega^n)\min(|\omega^n|,\omega^{max}).
$$
The relaxation factor in the first time step is equal to $\omega^{max}$.

This method improves convergence speed drastically compared to Gauss-Seidel iterations,
but even faster convergence can be obtained using quasi-Newton methods, which can be interpreted as using different relaxation factors for differnt Fourier modes of the ouput of the second solver.

### Settings

Besides the parameters required in the class `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`omega_max`|double|Maximal relaxation factor.

## IQNI

This coupled solver inherits from the class `GaussSeidel`.
The `type` for this coupled solver is `coupled_solvers.iqni`.

### Algorithm

The abbreviation IQNI refers to _interface quasi-Newton with inverse Jacobian_.
In this type of coupling iteration, the combination of the two solvers is seen as one system.
The input of the first solver $\mathcal{F}$ in iteration $k$ is denoted by $x^k$.
The output of this solver is transferred unchanged to the second solver.
The output of the second solver is denoted $\tilde{x}^k$.
The difference between input and output is called the residual $r^k=\tilde{x}^k-x^k$.

A residual operator $\mathcal{R}(x)$ is defined which return the residual $r^k$ as a function of $x^k$.
The goal is to find $x$ for which $\mathcal{R}(x)=0$, i.e. the root.
This system of non-linear equations is solved using Newton-Raphson iterations as follows
$$
\Delta x^k=\mathcal{R}'^{-1}(x^k)\Delta r,
$$
where $\mathcal{R}'$ is the Jacobian of $\mathcal{R}$ with respect to $x$, $\Delta x^k=x^{k+1}-x^k$ is the difference between the input of two subsequent iterations 
and $\Delta r^k=0-r^k=-r^k$ is the difference between the desired and the current residual
The iteration update can also be written as
$$
x^{k+1}=x^k-\mathcal{R}'^{-1}(x^k)r^k.
$$

However, this Jacobian is not accessible and therefore has to be approximated.
Instead of approximating $\mathcal{R}'$, solving the linear system can be avoided by approximating its inverse directly.

The approximation procedure will typically result in a low-rank Jacobian.
Whereas, a full rank Jacobian is required for the Newton-Raphson update to function properly.
Therefore, the inverse Jacobian of $\tilde{\mathcal{R}}$ is approximated.
This altered residual operator is defined as follows
$$
r^{k+1}=\tilde{\mathcal{R}}(\tilde{x}^{k+1})=\mathcal{R}(\tilde{x}^{k+1}-r^{k+1}).
$$
The inverse of both Jacobians are linked by
$$
\tilde{\mathcal{R}}'^{-1}=\mathcal{R}'^{-1}+I,
$$
where $\tilde{\mathcal{R}}'$ is the Jacobian of $\tilde{\mathcal{R}}$ with respect to $\tilde{x}$ and $I$ is the identity matrix.
For this Jacobian the following is valid
$$
\Delta \tilde{x}^k=\tilde{\mathcal{R}}'^{-1}(x^k)\Delta r^k,
$$
where $\Delta\tilde{x}^k=\tilde{x}^{k+1}-\tilde{x}^k$ is the difference between the the output of two subsequent iterations.
This Jacobian is also not known, but is approximated using a `model` and denoted by $\tilde{N}^k$.
The type of `model` and its settings are specified in the `settings` dictionary.
This model returns an estimation of $\Delta\tilde{x}^k$ given a value $\Delta r^k=-r^k$
$$
\Delta\tilde{x}^k=\tilde{N}^k \Delta r^k.
$$
Finally resulting in the update formula
$$
x^{k+1}=x^k+(\tilde{N}^k-I)\Delta r^k=x^k-\tilde{N}^k r^k+r^k.
$$

A symbolic schematic is given in the following figure.
![iqni](images/iteration_iqni.png "Residual operator iterations")

For more information with respect to the approximation of the Jacobian, refer to the [models documentation](models/models.md).
More information about residual operator methods can be found in [[1](#1)].

### Settings

Besides the parameters required in the class `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`model`|dict|Model component.
`omega`|double|Relaxation factor.

## IBQN

This coupled solver inherits from the class `GaussSeidel`.
The `type` for this coupled solver is `coupled solvers.ibqn`.

### Algorithm

The abbreviation IBQN refers to _interface block quasi-Newton_.
In type of coupling iteration, the system 
$$
\begin{cases}
    \mathcal{F}(x)-y=0
    \newline
    \mathcal{S}(y)-x=0
\end{cases}
$$
is solved in block form.
Again, $\mathcal{F}$ is the first solver with input $x$ and output $\tilde{y}$ 
and $\mathcal{S}$ is the second solver with input $y$ and output $\tilde{x}$
In this iteration scheme, the output of each solver is altered before being transferred to the other one.
Solving the system in block Newton-Raphson iterations results in
$$
	\begin{bmatrix}
		\mathcal{F}'(x) & -I
		\newline
		-I &   \mathcal{S}'(y)
	\end{bmatrix}
	\begin{bmatrix}
		\Delta x
		\newline
		\Delta y
	\end{bmatrix}
	=-
	\begin{bmatrix}
		\mathcal{F}(x)-y
		\newline
		\mathcal{S}(y)-x
	\end{bmatrix},
$$
where $\mathcal{F}'$ and $\mathcal{S}'$ denote the Jacobians of the first and second solver, respectively.
These Jacobians are, however, not accessible and are approximated using a `model` as specified in the `settings` dictionary.
To the fist and second solver correspond `model_f`, denoted here by $M_f$, and `model_s`, denoted bye $M_s$, respectively.
For example, `model_f` returns an estimation of $\Delta\tilde{y}^k=\tilde{y}^{k+1}-\tilde{y}^k$ given $\Delta x^k=x^{k+1}-x^k$
$$
\Delta\tilde{y}^k=M_f^k \Delta x^k.
$$

Solving for $x^{k+1}=x^k+\Delta x^k$ requires solving the system
$$
\left(I-M_s^k M_f^k\right)\Delta x^k
=\tilde{x}^k-x^k+M_s^k(\tilde{y}^k-y^k)
$$
for $\Delta x^k$.
This done matrix-free using the Generalized minimal residual method (GMRES).
Analogously, the input $y^{k+1}=y^k+\Delta y^k$ for the structural solver by solving 
$$
\left(I-M_f^{k+1}M_s^k\right)\Delta y^k
=\tilde{y}^{k+1}-y^k+M_f^{k+1}(\tilde{x}^k-x^{k+1})
$$
for $\Delta y^k$.

A symbolic schematic is given in the following figure.
![ibqn](images/iteration_ibqn.png "Block iterations")

The actual approximation of the Jacobian occurs with the same [models](models/models.md) as before.
More information about block methods can be found in [[1](#1)].

### Settings

Besides the parameters required in the class `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`absolute_tolerance_gmres`|double|Absolute tolerance used in the GMRES method.
`model_f`|dict|Model component corresponding to the first solver wrapper.
`model_s`|dict|Model component corresponding to the second solver wrapper.
`omega`|double|Relaxation factor.
`relative_tolerance_gmres`|double|Relative tolerance used in the GMRES method.

## Test single solver

The solver `test_single_solver` can be used to test new cases and solver settings.
The idea behind this component is to only test one of the two solvers, while the other one is replaced by a dummy.
This test environment inherits from the class `GaussSeidel`. 
The `type` for this coupled solver is `coupled_solvers.test_single_solver`.

### Dummy solver

To test only one solver, a dummy solver must be used.
Such a dummy solver is implemented by a test class in the file `dummy_solver.py`, which has to be on the same folder level as `run_simulation.py`.
Upon run-time an instance of this class is made.
The test class requires methods of the form `calculate_<variable>(x,y,z,n)`, with `<variable>` being the variable(s) required by the tested solver, e.g. `displacement`, `pressure` or `traction`.
How these variables are defined inside these methods, is up to the user.
However, the methods need to return the right format: a 3-element list or numpy array for vector variables and a 1-element list or numpy array for scalar variables.
Some examples are given in the example [test_single_solver](../../examples/test_single_solver/test_single_solver.md)
The test class name is provided in the JSON settings as a string.
If no test class is provided or the value `None` is used, zero input will be used.

### Settings

The JSON file requirements for the class `CoupledSolverTestSingleSolver` are different from the other coupled solvers in the sense that they only require 
the `type`, which is `coupled_solvers.test_single_solver`, the dictionary `test_settings` and the list `solver_wrappers` containing at least one `solver_wrapper`.
The keys for the `test_settings` dictionary are listed in alphabetical order below.

parameter|type|description
---:|:---:|---
`delta_t`|double|(optional) Time step size to be used in the test. Is optional as long as this value is defined in the `settings` dictionary. If a different value is defined in both dictionaries, the one defined in `test_settings` is chosen.
`name`|string|(optional) Name of the case used to store a pickle file with results. The pickle file will have the name `<name>_<test_solver_working_directory>.pickle`. If not provided, the value from `settings` is used or if `settings` is not present: "results".
`save_results`|bool|(optional) If `true` a [pickle](https://docs.python.org/3.9/library/pickle.html) file is stored containing some main [results](#save-results) as in `gauss_seidel`. If not provided, the value from `settings` is used or if `settings` is not present: `false`.
`solver_index`|int|Has a value 0 or 1 and indicates the solver that one wants to test. 0 indicates the first `solver_wrapper` that appears in the JSON-file, 1 the second one.
`test_class`|string|(optional) Refers to the class to use in the `dummy_solver.py`. If not provided or `None`, zero input will be used.
`timestep_start`|int|(optional) Time step to start from. If not provided the value defined in the `settings` dictionary is used. If the `settings` dictionary is not present, zero is used.

Other dictionaries, used for the actual calculation can be kept, but will not be used, with the possible exception of the `settings` dictionary.
The `settings` dictionary is used to look up `delta_t`, `timestep_start`, `save_results` and `name` if not provided in `test_settings`.
Note that `test_settings` has priority over the parameters defined in `settings`.
This means a calculation can be tested, by only adding the `test_settings` dictionary and changing the coupled solver `type` to `coupled_solvers.test_single_solver` and without altering anything else.
An illustration can be found in the example [test_single_solver](../../examples/test_single_solver/test_single_solver.md).

The working directory of the solver is copied to a new directory with the same name and a suffix `_testX` with `X` an integer starting from 0.  
As such, previous test solver working directories are not overwritten.
The optional pickle file, which saves some [results](#save-results), uses the name as specified by the JSON settings followed by an underscore and the solver test working directory.
As such, the pickle file can always be linked to the corresponding test working directory.

During run time, the norm of $x$ and $y$ are printed.
A residual does not exist here.
The arrays $x$ and $y$ do not have a physical meaning, but are the in- and output of the first solver,
which is typically the flow solver. Then, the vector $y$ will contain pressure and traction components for all points.
Nonetheless, these values can be useful to verify that the `solver_wrapper` runs.

The test environment `test_single_solver` tests only the `solver_wrapper` itself, no mapping is included.

## Save results
In each coupled solver, the `save_results` key can be set on `true`, in order to save some results into a [pickle](https://docs.python.org/3.9/library/pickle.html) file.
The key `name` dictates the name of this file as explained above.
The pickle file is saved after every time step and may be used by the postprocessing files included with the examples.
It contains a dictionary with the following keys (np-array refers to numpy array):

key|value type|description
---:|:---:|---
`solution_x`|numpy array|Contains the values of the vector $x$ (typically displacement) for every time step as columns of a 2 dimensional np-array. The initial value is included as well, such that the number of column will be one higher than the number of calculated time steps. The vector $x$ refers to the input of the first solver wrapper.
`solution_y`|numpy array|Contains the values of the vector $y$ (typically pressure and traction) for every time step as columns of a 2 dimensional np-array, similar to `solution_x`. The vector $y$ refers to the input of the second solver wrapper.
`interface_x`|interface|Interface object used as input for the first solver wrapper.
`interface_y`|interface|Interface object used as input for the second solver wrapper.
`iterations`|list|Contains the performed number of coupling iterations for every time step.
`time`|float|Equals the computation time between initialization and finalization of the final time step.
`residual`|list|Nested list, which contains for each time step a list of the residual after every iteration.
`delta_t`|float|Equals the used time step size.
`timestep_start`|int|Equals the used start time step.

In simulations with a large number of points on the interface and a very large number of time steps, this file might take up a larger amount of storage.

Finally, there is also a debug option, which can be activated by setting the boolean `self.debug` on `True` in the code file itself.
Then, the above information is stored every iteration. Additionaly, the residual vector is saved as well using the key `solution_r`, anologously to `solution_x` and `solution_y`.
Note, however, that the pickle file can be become very large in that case.

## References

<a id="1">[1]</a> 
[Delaissé N., Demeester T., Fauconnier D. and Degroote J., "Comparison of different quasi-Newton techniques for coupling of black box solvers", in ECCOMAS 2020, Proceedings, Paris, France, 2021.](http://hdl.handle.net/1854/LU-8685199)