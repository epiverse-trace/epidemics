# Create the connectivity matrix between all groups of `populations` using a gravity model

Create the connectivity matrix between all groups of `populations` using
a gravity model

## Usage

``` r
gravity_contact(populations, connectivity_matrix)
```

## Arguments

- populations:

  A list of `<population>` objects

- connectivity_matrix:

  A numeric matrix for the contact matrix between the elements of the
  populations

## Value

A numeric matrix combining the contact matrices in `<population>` and
the connectivity matrix using a gravity model
