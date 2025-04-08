## Improve `get_gene_expression_level` efficiency by performing adaptive sampling

`lik[k]` is either unimodal or uniform/degenerate (at least over large regions)
thus most iterations of the loop are unnecessary. We can use adaptive sampling to
find a region around the mode of `lik[k]`, using golden-section search eg, and
then sample more densely around that region. 

## Implement diagnostics for `Sanity()` fit

Check if the likelihood is unimodal, and if not, return a warning. The user can
then rerun it by adapting the range of the search.

## Generalize `calculateSanityDistance` to a generic method?
