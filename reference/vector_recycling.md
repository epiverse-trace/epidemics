# Vector recyclability checks and vector recycling

Internal functions to check whether vectors can be recycled and to
recycle vectors.

## Usage

``` r
.test_recyclable(x)

.recycle_vectors(x)
```

## Arguments

- x:

  A list of vectors to be checked for compliance with Tidyverse
  recycling rules, or to be recycled. When `x` is a list to be recycled,
  two cases are envisaged: all elements of `x` are expected to be of the
  same length, or, all but one element of `x` are scalars.

## Value

`.test_recyclable` returns a single logical value. `.recycle_vectors`
returns a list of vectors recycled to the same length following
Tidyverse recycling rules.

## Details

Note that `.test_recyclable()` will return vectors of unequal lengths if
`x` does not comply with length rules. This compliance is not enforced
as this function is only expected to be used internally, after a call to
`.test_recyclable()` (e.g. in a
[`stopifnot()`](https://rdrr.io/r/base/stopifnot.html)).

## Author

Tim Taylor
