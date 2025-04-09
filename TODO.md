## Improve `get_gene_expression_level` efficiency 

### Perform Adaptive Sampling

`lik[k]` is either unimodal or uniform/degenerate (at least over large regions)
thus most iterations of the loop are unnecessary. We can use adaptive sampling to
find a region around the mode of `lik[k]`, using golden-section search eg, and
then sample more densely around that region. 

### Use Newton-Raphson instead of bisection

Most or many of the counts are 0, most time is spent solving the `get_epsilon_2`
We can improve the efficiency by using an Newton-Raphson solver with:

```r
dL <- function(sigma, d, frac) sigma * (d + 0.5 * sigma) / v + gene_size * frac * expm1(sigma) - 0.5
dL_deriv <- function(sigma, d, frac) (d + sigma) / v + gene_size * frac * exp(sigma)
```

## Generalize `calculateSanityDistance` to a generic method?
