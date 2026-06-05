# Construct a new population for an epidemic model

Construct a new population for an epidemic model

Check whether an object is a `<population>`

## Usage

``` r
population(
  name = NA_character_,
  contact_matrix,
  demography_vector,
  initial_conditions
)

is_population(x)
```

## Arguments

- name:

  Optional string for the population name.

- contact_matrix:

  A matrix giving the contacts between the demographic groups in the
  population. Must be a square matrix.

- demography_vector:

  A vector of the sizes of each demographic group. Must have the same
  length as the dimensions of the contact matrix.

- initial_conditions:

  Matrix representing the initial proportions of each demographic group
  in the four model compartments: 'susceptible', 'infected/infectious',
  'recovered', and 'vaccinated'. Must have as many rows as the number of
  demographic groups. Each compartment is represented in the columns of
  the matrix, so that the element \\M\_{ij}\\ represents the proportion
  of individuals of demographic group \\i\\ in compartment \\j\\ .

- x:

  An object to be checked as a valid population.

## Value

An object of the `<population>` S3 class.

`is_population()` returns a logical for whether the object is a
`<population>`.

## Examples

``` r
uk_pop <- population(
  name = "UK population",
  contact_matrix = matrix(1),
  demography_vector = 67e6,
  initial_conditions = matrix(
    c(0.9999, 0.0001, 0, 0),
    nrow = 1, ncol = 4
  )
)

# print to check
uk_pop
#> <population> object
#> 
#>  Population name: 
#> "UK population"
#> 
#>  Demography 
#> Dem. grp. 1: 67,000,000 (100%)
#> 
#>  Contact matrix 
#>              Dem. grp. 1:
#> Dem. grp. 1:            1
#> 
#>  Initial Conditions 
#>        [,1]  [,2] [,3] [,4]
#> [1,] 0.9999 1e-04    0    0

# check for class <population>
is_population(uk_pop)
#> [1] TRUE
```
