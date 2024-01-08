# Convergence criteria

Convergence criteria are an essential part of numerical tools. They should be chosen wisely in order to obtain a reasonably
accurate solution without performing more iterations than needed. This documentation describes how the user can practically
assemble a set of convergence criteria.

## Types

### Iteration limit

The `type` `convergence_criterion.iteration_limit` is satisfied if the number of coupling iterations equals or is larger than a predefined maximum.

The `settings` dictionary contains one entry:

| parameter | type | description                   |
|----------:|:----:|-------------------------------|
| `maximum` | int  | Maximum number of iterations. |

### Absolute norm

The `type` `convergence_criterion.absolute_norm` is satisfied if the $p$-norm of the residual in the last coupling iteration
is smaller than a predefined `tolerance`. 
More information on how the residual is calculated can be found in the [coupled solvers documentation](../coupled_solvers/coupled_solvers.md).
The $p$-norm of the residual $r$ is defined as

$$
\Vert r \Vert_p = \left(\sum_{i=1}^n \vert r_{i}\vert^p\right)^{1/p} ,
$$

where $r_i$ is the $i$-th component of the residual.

The `settings` are as follows:

|   parameter |  type  | description                  |
|------------:|:------:|------------------------------|
|     `order` |  int   | Order $p$ of the norm.       |
| `tolerance` | double | Limit value for convergence. |

### Relative norm

The `type` `convergence_criterion.relative_norm` is completely analogous to the absolute norm. Instead of the norm of 
the last residual being smaller than a set `tolerance`, now the ratio norm of the residual of the last coupling iteration to the 
norm of the residual of the first coupling iteration is compared to a `tolerance`. Zero divisions are avoided internally by
comparing the norm of the residual of the first coupling iteration to the machine limit for floating points, i.e. the smallest number
different from zero. In case this initial norm is too small, an exception will be raised. Should this happen, it is advised to opt
for the absolute norm criterion instead of relative norm.

These are the same as for `convergence_criterion.absolute_norm`.

### Solver convergence

The `type` `convergence_criterion.solver_convergence` is a convergence criterion that uses the convergence of a solver to determine the convergence of the coupling loop.
The criterion is considered fulfilled once the solver converges in its first solver iteration.
Not all solvers allow for this feature, see their respective documentation.

Typically, a `solver_convergence` criterion is used for each solver, combined with an `and` criterion (see [Combining multiple convergence criteria](#combining-multiple-convergence-criteria)).
This way, the coupling loop is considered converged when all solvers converge in their first solver iteration [[1](#1)].

The `settings` dictionary contains one entry:

|      parameter | type | description                                                                                                                 |
|---------------:|:----:|-----------------------------------------------------------------------------------------------------------------------------|
| `solver_index` | int  | Index of the solverwrapper in the `solver_wrappers` list of the coupled solver to which this convergence criterion belongs. |

## Combining multiple convergence criteria

In most cases, it is wise to combine two criteria: one to ensure a high enough accuracy and an iteration limit in order
to break loops that are not converging fast enough. In that case, the criteria are combined via `or` or `and` statements.
In that case, the `type` is set to `convergence_criterion.or` (alternatively, `convergence_criterion.and`) and the `settings`
contain a `criteria_list` that contains single criteria in the same way as described above. 

In the following example, the `iteration_limit` and `relative_norm` criteria are combined using an `or` statement.
Note that the number of criteria is not limited to two.
Moreover, the `or` and `and` statements can be combined multiple times, if needed.

```json
{
  "type": "convergence_criteria.or",
  "settings": {
    "criteria_list": [
      {
        "type": "convergence_criteria.iteration_limit",
        "settings": {
          "maximum": 20
        }
      },
      {
        "type": "convergence_criteria.relative_norm",
        "settings": {
          "order": 2,
          "tolerance": 1e-3
        }
      }
    ]
  }
}
```

An example of the use of the `and` criterion is with use of the `solver_convergence` criterion.

```json
{
  "type": "convergence_criteria.and",
  "settings": {
    "criteria_list": [
      {
        "type": "convergence_criteria.solver_convergence",
        "settings": {
          "solver_index": 0
        }
      },
      {
        "type": "convergence_criteria.solver_convergence",
        "settings": {
          "solver_index": 1
        }
      }
    ]
  }
}
```

## Dummy convergence criterion

This dummy convergence criterion can be used in the [explicit](../coupled_solvers.md#explicit) or [one-way](../coupled_solvers.md#one-way) coupled solvers, which don't require a convergence criterion.

If the use of a dummy convergence criterion is allowed, the `type` (`convergence_criteria.dummy_convergence_criterion`) can be written explicitly or omitted. No `settings` are required.
Note that the key `convergence_criterion` is still required.

## References

<a id="1">[1]</a> 
[Spenke T., Delaissé N., Degroote J. and Hosters, N., "On the number of subproblem iterations per coupling step in partitioned fluid-structure interaction simulations", Preprint on ArXiv, 2023.](
https://doi.org/10.48550/arXiv.2303.08513)
