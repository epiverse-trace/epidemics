# Create the connectivity matrix between all groups of `populations`

Create the connectivity matrix between all groups of `populations`

## Usage

``` r
combine_contact(populations, connectivity_matrix)
```

## Arguments

- populations:

  A list of `<population>` objects

- connectivity_matrix:

  A numeric matrix for the proportion of contacts between the elements
  of the populations

## Value

A numeric matrix combining the contact matrices in `<population>` and
the connectivity matrix.
