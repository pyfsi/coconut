# CoupledSolvers

This documentation describes the coupled solvers which are available.
A coupled solver refers to a coupling algortihm used to couple the two solvers.
Some coupled solvers implement one or more models.
All coupled solvers inherit from `gauss_seidel`.

In the parameter JSON file, the dictionary `coupled_solver` holds the dictionaries `type` and `settings`,
but also `predictor`, `convergence_criterion` and `solver_wrappers`.

## GaussSeidel

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
`time_step_start`|int|Time step number to (re)start a transient FSI calculation. If 0 is given, the simulation starts from scratch. Otherwise, the code looks for the relevant files to start from the corresponding time step. Not every solver implements restart, see the solver documentation for more information. For a steady simulation this should be 0.

`timestep_start` and `delta_t` are necessary parameters (also in a steady simulation), but can also defined in the solverwrapper directly (e.g. for standalone testing).
If they are defined both here and in the solver wrapper, then the former value is used and a warning is printed.

**These parameters should also be specified for the coupled solvers inheriting from `GaussSeidel`.**

## Relaxation

This coupled solvers inherits from `GaussSeidel`.

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

This coupled solvers inherits from `GaussSeidel`.

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

This coupled solvers inherits from `GaussSeidel`.

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

For more information refer to .

### Settings

Besides the parameters required in `GaussSeidel`, the following parameters need to be included in the `settings` dictionary.
They are listed in alphabetical order.

parameter|type|description
---:|:---:|---
`model`|dict|Model component.
`omega`|double|Relaxation factor.

## IBQN

This coupled solvers inherits from `GaussSeidel`.

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
