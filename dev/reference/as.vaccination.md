# Convert a list to a vaccination object

Convert a list to a vaccination object

## Usage

``` r
as.vaccination(x)
```

## Arguments

- x:

  A list, or an object that inherits from a list.

## Value

A
[vaccination](https://epiverse-trace.github.io/epidemics/dev/reference/vaccination.md)
class object.

## Examples

``` r
# prepare a list
vax <- list(
  name = "vax_regime",
  time_begin = matrix(1),
  time_end = matrix(100),
  nu = matrix(0.001)
)

as.vaccination(vax)
#> <vaccination> object
#> 
#>  Vaccination name: 
#> "vax_regime"
#> 
#>  Begins at: 
#>      dose_1
#> [1,]      1
#> 
#>  Ends at: 
#>      dose_1
#> [1,]    100
#> 
#>  Vaccination rate: 
#>      dose_1
#> [1,]  0.001
```
