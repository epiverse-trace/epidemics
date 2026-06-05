# Constructor for the super-class and sub-classes

Constructor for the super-class and sub-classes

Constructor for a new \<contacts_intervention\>

Constructor for a new \<rate_intervention\>

## Usage

``` r
new_intervention(
  name = NA_character_,
  time_begin,
  time_end,
  reduction,
  ...,
  class
)

new_contacts_intervention(name, time_begin, time_end, reduction)

new_rate_intervention(name, time_begin, time_end, reduction)
```

## Arguments

- name:

  String for the name of the intervention.

- time_begin:

  A matrix with a single element giving the start time of the
  intervention.

- time_end:

  A matrix with a single element giving the end time of the
  intervention.

- reduction:

  A matrix with a single column and as many rows as there are
  demographic groups to be targeted by the intervention. See
  [`intervention()`](https://epiverse-trace.github.io/epidemics/dev/reference/intervention.md)
  for details on the types of intervention that can be created and the
  requirements for the type of `reduction`.

- ...:

  Any other parameters to be passed to the constructor.

- class:

  A string giving the type of the intervention; used to generate
  intervention sub-classes.

## Value

`new_intervention()` returns an object that inherits from the
`<intervention>` class. `new_contacts_intervention()` returns an object
that is a sub-class of `<intervention>` called
`<contacts_intervention>`. `new_rate_intervention()` returns an object
that is a sub-class of `<intervention>` called `<rate_intervention>`.
