# Convergence criteria

Convergence criteria are an essential part of numerical tools. They should be chosen wisely in order to obtain a reasonably
accurate solution without performing more iterations than needed. This documentation describes how the user can practically
assemble a set of convergence criteria.

## Types

### Iteration limit

The `type` `convergence_criterion.iteration_limit` is satisfied if the number of coupling iterations equals or is larger than a predefined maximum.

#### Settings

The `settings` dictionary contains one entry:

parameter|type|description
---:|:---:|---
`maximum`|int|Maximum number of iterations.

### Absolute norm

The `type` `convergence_criterion.absolute_norm` is satisfied if the $p$-norm of the residual in the last coupling iteration
is smaller than a predefined `tolerance`. The $p$-norm of the residual $r$ is defined as

$$
\Vert r \Vert_p = \left(\sum_{i=1}^n \vert r_{i}\vert^p\right)^{1/p} ,
$$

where $r_i$ is the $i$-th component of the residual.

#### Settings

The `settings` are as follows:

parameter|type|description
---:|:---:|---
`tolerance`|double|Limit value for convergence.
`order`|int|Order $p$ of the norm.

### Relative norm

This `type` (`convergence_criterion.relative_norm`) is completely analogous to the absolute norm. Instead of the norm of 
the last residual being smaller than a set `tolerance`, now the ratio norm of the residual of the last coupling iteration to the 
norm of the residual of the first coupling iteration is compared to a `tolerance`. Zero divisions are avoided internally by
comparing the norm of the residual of the first coupling iteration to the machine limit for floating points, i.e. the smallest number
different from zero. In case this initial norm is too small, an exception will be raised. Should this happen, it is advised to opt
for the absolute norm criterion instead of relative norm.

#### Settings

These are the same as for `convergence_criterion.absolute_norm`.

## Combining multiple convergence criteria

In most cases, it is wise to combine two criteria: one to ensure a high enough accuracy and an iteration limit in order
to break loops that are not converging fast enough. In that case, the criteria are combined via `or` or `and` statements.
In that case, the `type` is set to `convergence_criterion.or` (alternatively, `convergence_criterion.and`) and the `settings`
contain a `criteria_list` that contains single criteria in the same way as described above. 

In the following example, the `iteration_limit` and `relative_norm` criteria are combined using an `or` statement. 

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
          "tolerance": 1e-3,
          "order": 2
        }
      }
    ]
  }
}
```


