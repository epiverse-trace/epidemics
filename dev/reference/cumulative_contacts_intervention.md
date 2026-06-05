# Calculate the Cumulative Effect of Interventions on Social Contacts

Calculate the Cumulative Effect of Interventions on Social Contacts

Scale a Contact Matrix by all Active Interventions on Social Contacts

## Usage

``` r
cumulative_contacts_intervention(t, time_begin, time_end, reduction)

intervention_on_cm(t, cm, time_begin, time_end, cr)
```

## Arguments

- t:

  The current time.

- time_begin:

  A numeric vector of the start times of all interventions being
  modelled.

- time_end:

  A numeric vector of the end times of all interventions being modelled.
  Must be the same length as `time_begin`.

- reduction:

  A numeric matrix where rows give the effect of interventions on each
  demographic group, and the columns give the proportion reduction in
  contacts. When two interventions overlap, the proportions are *added*,
  for a maximum possible value of 1.0 (i.e., no contacts).

- cm:

  A numeric matrix of social contacts between demographic groups.

## Value

`cumulative_contacts_intervention()` returns a numeric vector of the
proportion reduction in contacts for each demographic group.

`intervention_on_cm()` returns the contact matrix `cm` scaled by the
cumulative effect of any active interventions.
