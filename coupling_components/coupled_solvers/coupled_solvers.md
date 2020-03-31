# CoupledSolvers

This documentation describes the coupled solvers which are available.
A coupled solver refers to a coupling algortihm used to couple the two solvers.
Some coupled solvers implement one or more models.
All coupled solvers inherit from `gauss_seidel`.

In the parameter JSON file, the dictionnary `coupled_solver` holds the dictionnaries `type` and `settings`,
but also `predictor`, `convergence_criterion` and `solver_wrappers`.

## Gauss-Seidel

The following parameters, listed in alphabetical order, need to be specified in the `solver_parameters.json` file in the `working_directoyr`.
Care should be taken that the values of `d`, `e`, `h`, `l` and `rhof` match the corresponding values of the flow solver.

parameter|type|description
---:|:---:|---
`axial_offset`|double|(optional) Distance over which tube is displaced axially in the coordinate system. If not provided, the value of this parameter is 0.
`d`|double|Nominal diameter of the tube.
`e`|double|Modulus of elasticity of the tube wall.

## Relaxation

## Aitken

## IQN-I

## IBQN
