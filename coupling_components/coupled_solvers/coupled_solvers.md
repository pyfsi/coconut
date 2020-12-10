# Coupled solvers

This documentation describes the coupled solvers that are available.
A coupled solver refers to a coupling algortihm used to couple the two solvers.
Some coupled solvers implement one or more models.
An odd one out is `test_single_solver` which allows to test only one solver by coupling it to a dummy solver. 
All coupled solvers inherit from `gauss_seidel`.

In the parameter JSON file, the dictionary `coupled_solver` holds the `type` and the dictionary `settings`,
but also the dictionary `predictor`, the dictionary `convergence_criterion` and the list `solver_wrappers` containing 2 dictionaries, one for each `solver_wrapper`.
More information on these last three can be found in respectively [predictors](../predictors/predictors.md),
[convergence_criteria](../convergence_criteria/convergence_criteria.md) and the `solver_wrappers` documentation.


## GaussSeidel

The `type` for this `coupled_solver` is `coupled_solvers.gauss_seidel`.

### Algorithm

Gauss-Seidel iterations are the simplest way of coupling two solvers.
In this coupling method fixed point iterations are performed as the output of one solver is given to the other one without adjustment.

Basic schematic is given in the following figure.
![gauss_seidel](images/iteration_gauss_seidel.png "Gauss-Seidel iterations")

Here, $\mathcal{F}$ is the first solver with input $x$ and output $\tilde{y}$ 
and $\mathcal{S}$ is the second solver with input $y$ and output $\tilde{x}$

Gauss-Seidel iterations are very simple, but unstable for cases with high added mass.

### Settings

The following parameters need to be included in the `settings` dictionary.
Here they are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`delta_t`|double|Fixed time step size used in both solvers. For a steady simulation typically a value of 1 is taken.
`name`|string|(optional) Name of the case used to store a pickle file with results. If not provided 'results' is used.
`save_results`|bool|(optional) Default: false. If true a pickle file is stored containing some main results: the value of the displacement and load on the interface for every time step, interface objects used by the two solvers, the number of coupling iterations per time step, the total elapsed time for the calculation, the residual after every coupling iteration and the values of `delta_t` and `timestep_start`.<br> This file is used by the postprocessing files included with the examples.
`time_step_start`|int|Time step number to (re)start a transient FSI calculation. If 0 is given, the simulation starts from scratch. Otherwise, the code looks for the relevant files to start from the corresponding time step. Not every solver implements restart, see the solver documentation for more information. For a steady simulation this should be 0.

`timestep_start` and `delta_t` are necessary parameters (also in a steady simulation), but can also defined in the solverwrapper directly (e.g. for standalone testing).
If they are defined both here and in the solver wrapper, then the former value is used and a warning is printed.

**These parameters should also be specified for the coupled solvers inheriting from `GaussSeidel`.**

## Relaxation

This `coupled_solver` inherits from `GaussSeidel`.
The `type` for this `coupled_solver` is `coupled_solvers.relaxation`.

### Algorithm

Gauss-Seidel iterations are very simple, but are unstable for cases with high added mass.
An approach to mitigate this is applying relaxation, also called simple mixing.
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

Here, $\mathcal{F}$ is the first solver with input $x$ and output $\tilde{y}$ 
and $\mathcal{S}$ is the second solver with input $y$ and output $\tilde{x}$

A lower $\omega$ increases stability, but decreases convergence speed.

This method is again quite simple, but still often unstable for cases with high added mass.

### Settings

Besides the parameters required in `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`omega`|double|Relaxation factor.

## Aitken

This `coupled_solver` inherits from `GaussSeidel`.
The `type` for this `coupled_solver` is `coupled_solvers.aitken`.

### Algorithm

Gauss-Seidel iterations are very simple, but are unstable for cases with high added mass.
An approach to mitigate this is applying relaxation.
However, the choice of a relaxation factor can be difficult.
In this coupling method a dynamic relaxation factor is used.
The output of first solver is still given to the other one without adjustment, 
but the output of the second solver is relaxed as follows:
$$
x^{k+1}=(1-\omega^k)x^k+\omega^k\tilde{x}^k=x^k+\omega^k r^k
$$
with $x^k$ and $\tilde{x}^k$, respectively the input for the first solver and the output of the second solver in iteration $k$.
The difference between both is called the residual $r^k=\tilde{x}^k-x^k$.
The mixing or relaxation factor is $\omega^k$ is dynamic in the sense that its value changes between iterations.

A symbolic schematic is given in the following figure.
![aitken](images/iteration_aitken.png "Aitken iterations")

Here, $\mathcal{F}$ is the first solver with input $x$ and output $\tilde{y}$ 
and $\mathcal{S}$ is the second solver with input $y$ and output $\tilde{x}$

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
but even faster convergence can be obtained using quasi-Newton methods which take into account sensitivities.

### Settings

Besides the parameters required in `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`omega_max`|double|Maximal relaxation factor.

## IQNI

This `coupled_solver` inherits from `GaussSeidel`.
The `type` for this `coupled_solver` is `coupled_solvers.iqni`.

### Algorithm

The abbreviation IQNI refers to _interface quasi-Newton with inverse Jacobian_.
In this type of coupling iteration, the combination of the two solvers is seen as one system.
The input of the first solver in iteration $k$ is denoted by $x^k$.
The output of this solver is transferred unchanged to the second solver.
The output of the second solver is denoted $\tilde{x}^k$.
The difference between input and output is called the residual $r^k=\tilde{x}^k-x^k$.

A residual operator $\mathcal{R}(x)$ is defined which return the residual $r^k$ as a function of $x^k$.
The goal is to find $x$ for which $\mathcal{R}(x)=0$.
This system of non-linear equations is solved using Newton-Raphson iterations as follows
$$
x^{k+1}=x^k-\mathcal{R}'^{-1}(x^k)r^k,
$$
where $\mathcal{R}'$ is the Jacobian of $\mathcal{R}$. 
However, this Jacobian is not accessible and is therefore approximated.

Note that the iteration update can also be written as
$$
\Delta x^k=\mathcal{R}'^{-1}(x^k)\Delta r^k;
$$
where $\Delta x^k=x^{k+1}-x^k$ is the difference between the input of two subsequent iterations 
and $\Delta r^k=0-r^k=-r^k$ is the difference between the desired and the current residual. 
Likewise, $\Delta\tilde{x}^k=\tilde{x}^{k+1}-\tilde{x}^k$ is the difference between the the output of two subsequent iterations.

Instead of approximating $\mathcal{R}'^{-1}$ directly, the inverse Jacobian of $\tilde{\mathcal{R}}$ is approximated.
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
\Delta \tilde{x}^k=\tilde{\mathcal{R}}'^{-1}(x^k)\Delta r^k.
$$
This Jacobian is also not known, but is approximated using a `model` denoted by $\tilde{N}^k$.
The type of `model` and its settings are specified in the `settings` dictionary.
This model returns an estimation of $\Delta\tilde{x}^k$ given $\Delta r^k=-r^k$
$$
\Delta\tilde{x}^k=\tilde{N}^k \Delta r^k.
$$
Finally resulting in the update formula
$$
x^{k+1}=x^k+(\tilde{N}^k-I)\Delta r^k=x^k-\tilde{N}^k r^k+r^k.
$$

A symbolic schematic is given in the following figure.
![iqni](images/iteration_iqni.png "Residual operator iterations")

Here, $\mathcal{F}$ is the first solver with input $x$ and output $\tilde{y}$ 
and $\mathcal{S}$ is the second solver with input $y$ and output $\tilde{x}$

For more information refer to [models](models/models.md).

### Settings

Besides the parameters required in `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`model`|dict|Model component.
`omega`|double|Relaxation factor.

## IBQN

This `coupled_solver` inherits from `GaussSeidel`.
The `type` for this `coupled_solver` is `coupled_solvers.ibqn`.

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
Here, $\mathcal{F}$ is the first solver with input $x$ and output $\tilde{y}$ 
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
These Jacobians are however not accessible and are approximated using a `model` as specified in the `settings` dictionary.
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
This done matrix-free using a GMRES method.
Analogously, the input $y^{k+1}=y^k+\Delta y^k$ for the structural solver by solving 
$$
\left(I-M_f^{k+1}M_s^k\right)\Delta y^k
=\tilde{y}^{k+1}-y^k+M_f^{k+1}(\tilde{x}^k-x^{k+1})
$$
for $\Delta y^k$.

A symbolic schematic is given in the following figure.
![ibqn](images/iteration_ibqn.png "Block iterations")

For more information refer to .

### Settings

Besides the parameters required in `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
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
This test environment inherits from `GaussSeidel`. 
The `type` for this `coupled_solver` is `coupled_solvers.test_single_solver`.

### Dummy solver

To test only one solver, a dummy solver must be used.
Such a dummy solver is implemented by a test class in the file `dummy_solver.py`, which has to be on the same folder level as `run_simulation.py`.
Upon run-time an instance of this class is made.
The test class requires methods of the form `calculate_<variable>(x,y,z,n)`, with `<variable>` being the variable(s) required by the tested solver, e.g. `DISPLACEMENT`, `PRESSURE` or `TRACTION`.
How these variables are defined inside these methods, is up to the user.
However, the methods need to return the right format: a 3-element list of floats for vector variables and a single float for scalar variables.
Some examples are given in the example [test_single_solver](../../examples/test_single_solver/test_single_solver.md)
The test class name is provided in the JSON settings as a string.
If no test class is provided or the value `None` is used, zero input will be used.

### Settings

The JSON settings for this test environment `test_single_solver` are different from the other `coupled_solvers` in the sense that they only require 
the `type`, which is `coupled_solvers.test_single_solver`, the dictionary `test_settings` and the list `solver_wrappers` containing at least one `solver_wrapper`.
The possibilities for the `test_settings` dictionary are listed in alphabetical order below.

parameter|type|description
---:|:---:|---
`delta_t`|double|(optional) Time step size to be used in the test. Is optional as long as this value is defined in the `settings` dictionary. If a different value is defined in both dictionaries, the one defined in `test_settings` is chosen.
`name`|string|(optional) Name of the case used to store a pickle file with results. If not provided 'results' is used. The pickle file will have the name `<name>_<test_solver_working_directory>.pickle`.
`save_results`|bool|(optional) Default: false. If true a pickle file is stored containing some main results as in `gauss_seidel`.
`solver_index`|int|Has a value 0 or 1 and indicates the solver that one wants to test. 0 indicates the first `solver_wrapper` that appears in the JSON-file, 1 the second one.
`test_class`|string|(optional) Refers to the class to use in the `dummy_solver.py`. If not provided or `None`, zero input will be used.
`timestep_start`|int|(optional) Time step to start from. If not provided the value defined in the `settings` dictionary is used. If the `settings` dictionary is not present, zero is used.

Other dictionaries, used for the actual calculation can be kept, but will not be used, with the possible exception of the `settings` dictionary.
The `settings` dictionary is used to look up `delta_t`, `timestep_start`, `save_results` and `name` if not provided in `test_settings`.
Note that `test_settings` has priority over the parameters defined in `settings`.
This means a calculation can be tested, by only adding the `test_settings` dictionary and changing the `coupled_solver` type to `coupled_solvers.test_single_solver` and without altering anything else.
An illustration can be found in the example [test_single_solver](../../examples/test_single_solver/test_single_solver.md).

The working directory of the solver is copied to a new directory with the same name and a suffix `_testX` with `X` an integer starting from 0.  
As such, previous test solver working directories are not overwritten.
The optional pickle file saving some results uses the name as specified by the JSON settings followed by an underscore and the solver test working directory.
As such the pickle file always belongs to the corresponding test working directory.

During run time, the norm of $x$ and $y$ are printed.
A residual does not exist here.
The arrays $x$ and $y$ do not have a physical meaning but are the in- and output of the solver:
the output of a typical flow solver, for example, will contain pressure and traction components for all points.
For the first solver in the `solver_wrappers` list, $x$ and $y$ are respectively in- and output,
whereas for the second, $y$ is the input an $x$ the output.
Nonetheless, these values can be useful to verify that the `solver_wrapper` runs.

The test environment `test_single_solver` tests only the `solver_wrapper` itself, no mapping is included.