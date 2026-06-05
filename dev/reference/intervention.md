# Create an intervention for an epidemic model

Prepare an object of the `<intervention>` super-class that specifies a
modification of the model parameters.

A `<contacts_intervention>` is used to simulate a non-pharmaceutical
intervention (NPI) regime that reduces the population's social contacts.

A `<rate_intervention>` is used to simulate a reduction in the model's
rate parameters (such as the transmission rate \\\beta\\), and can be
used to represent pharmaceutical interventions such as improved
treatment, but also NPIs such as wearing masks.

Interventions have a single start and end time that applies to all
demographic groups in the population, but can have groups-specific
effects on the reduction of contacts.

Combine `<intervention>`-inheriting objects to create sequential or
overlapping intervention regimes using
[`c()`](https://rdrr.io/r/base/c.html) on two or more
`<intervention>`-inheriting objects.

## Usage

``` r
intervention(name = NA_character_, type, time_begin, time_end, reduction)

is_intervention(x)

is_contacts_intervention(x)

is_rate_intervention(x)

# S3 method for class 'contacts_intervention'
c(x, ...)

# S3 method for class 'rate_intervention'
c(x, ...)
```

## Arguments

- name:

  String for the name of the intervention.

- type:

  String for the type of intervention. May be one of `"contacts"` or
  `"rate"`, for a `<contacts_intervention>` or `<rate_intervention>`
  respectively.

- time_begin:

  Single number for the start time of the intervention.

- time_end:

  Single number for the end time of the intervention.

- reduction:

  For `<contacts_intervention>`s, a matrix with as many rows as the
  number of demographic groups in the type population, and a single
  column. Each element gives the group-specific proportion reduction in
  contacts.

  For `<rate_intervention>`s, a single number giving the proportion
  reduction in a model parameter contacts.

  See details for how [`c()`](https://rdrr.io/r/base/c.html) can be used
  to combine interventions of the same sub-class.

- x:

  An `<intervention>` object, or an object to be checked as an
  `<intervention>` object.

- ...:

  intervention objects to combine with `x` to create a multi-dose
  `<intervention>` object.

## Value

An object of the `<intervention>` S3 super-class, with possible
sub-classes `<contact_intervention>` and `<rate_intervention>`.

Concatenating two or more `<intervention>`-inheriting objects using
[`c()`](https://rdrr.io/r/base/c.html) also returns a
`<intervention>`-inheriting object of the same sub-class. This object
holds the intervention-specific start and end times, and reductions
specified by all the constituent intervention actions (by demographic
group if an intervention on contacts).

The combined effect of these actions on the population is handled
internally by epidemic model functions.

A "null" intervention generated using
`.no_contacts_intervention(population)` or `.no_rate_intervention()`
returns a `<intervention>` of the appropriate sub-class that has its
start and end times, and its effect all set to 0.0.

`is_intervention()`, `is_contacts_intervention()`, and
`is_rate_intervention()` each return a logical value for whether the
object is of the `<intervention>`, `<contacts_intervention>`, or
`<rate_intervention>` class, respectively.

## Details

Epidemic models that can accommodate interventions on contacts are able
to accommodate any number of interventions with different start and end
times and different group-specific effects.

Epidemic models that can accommodate interventions on rates are also
able to accommodate any number of interventions with different start and
end times, but with only a uniform effect on the relevant rate.

When multiple contact interventions are combined using
[`c()`](https://rdrr.io/r/base/c.html), the reduction in contacts is
stacked column wise to form a matrix \\\[i, j\]\\.

When multiple rate interventions are combined using
[`c()`](https://rdrr.io/r/base/c.html), the reduction in the rate is
concatenated into a vector of the same length as the number of
interventions.

Models such as model_default() are set up to treat interventions with
overlapping periods (i.e., overlap between the time when they are active
) as having an *additive effect* on contact or rate reductions.

For contact reductions, the group-specific effect of \\J\\ overlapping
interventions is thus a vector \\\sum\_{j = 1}^J x\_{ij}\\, for each
demographic group \\i\\. This is handled internally by the epidemic
model code. For example, a contact reduction matrix for two perfectly
overlapping interventions (\\J = 2\\) with different effects across
three demographic groups (\\I = 3\\) would be represented as:
\\\begin{bmatrix}0.1 & 0.05\\0.1 & 0.1\\0.1 & 0.0\end{bmatrix}\\ In
epidemic models, the cumulative group-specific effect when both
interventions are active would be \\(0.15, 0.2, 0.1)\\.

For rate reductions, the effect of overlapping interventions that reduce
a particular rate is also considered to be additive.

## Examples

``` r
# assuming a population with two age groups, 0 -- 18, and 18+
# an example in which schools are closed for 30 days (or other time units)
close_schools <- intervention(
  name = "close schools",
  type = "contacts",
  time_begin = 50,
  time_end = 80,
  reduction = matrix(c(0.5, 0.01)) # reduces contacts differentially
)
close_schools
#> <contacts_intervention> object
#> 
#>  Intervention name: 
#> "close schools"
#> 
#>  Begins at: 
#> [1] 50
#> 
#>  Ends at: 
#> [1] 80
#> 
#>  Reduction: 
#>              Interv. 1
#> Demo. grp. 1      0.50
#> Demo. grp. 2      0.01

# Check for intervention class
is_contacts_intervention(close_schools)
#> [1] TRUE

# Concatenating interventions
# create first intervention
npi_1 <- intervention(
  type = "contacts",
  time_begin = 30,
  time_end = 60,
  reduction = matrix(0.1)
)

# second intervention
npi_2 <- intervention(
  type = "contacts",
  time_begin = 45,
  time_end = 75,
  reduction = matrix(0.1)
)

c(npi_1, npi_2)
#> <contacts_intervention> object
#> 
#>  Intervention name: 
#> NA
#> 
#>  Begins at: 
#>      npi_1 npi_2
#> [1,]    30    45
#> 
#>  Ends at: 
#>      npi_1 npi_2
#> [1,]    60    75
#> 
#>  Reduction: 
#>              Interv. 1 Interv. 2
#> Demo. grp. 1       0.1       0.1
```
