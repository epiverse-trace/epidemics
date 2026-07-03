# Convert a list to a intervention object

Convert a list to a intervention object

## Usage

``` r
as.intervention(x, type = c("contacts", "rate"))
```

## Arguments

- x:

  A list, or an object that inherits from a list.

- type:

  A string for the type of intervention: `"contacts"` for a
  `<contact_intervention>` or `"rate"` for a `<rate_intervention>`.

## Value

A
[intervention](https://epiverse-trace.github.io/epidemics/reference/intervention.md)
class object.

## Examples

``` r
# prepare a list
npi <- list(
  name = "npi",
  type = "contacts",
  time_begin = 30,
  time_end = 60,
  reduction = rep(0.1, 3)
)

as.intervention(npi)
#> <contacts_intervention> object
#> 
#>  Intervention name: 
#> "npi"
#> 
#>  Begins at: 
#> [1] 30
#> 
#>  Ends at: 
#> [1] 60
#> 
#>  Reduction: 
#>              Interv. 1
#> Demo. grp. 1       0.1
#> Demo. grp. 2       0.1
#> Demo. grp. 3       0.1
```
