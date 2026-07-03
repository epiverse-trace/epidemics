# Return ODE model output as a data.table

Return ODE model output as a data.table

Return Ebola model output as a table

## Usage

``` r
.output_to_df(output, population, compartments)

.output_to_df_ebola(output, population, compartments)
```

## Arguments

- output:

  The model output, which must be a two element list (for epidemic)
  models, with the names "x" and "time", where "x" represents the
  condition of each compartment at each timestep in "time".

- population:

  A `<population>` object corresponding to the population used in the
  epidemic model. The `<population>` object is used to generate the
  names of the demographic groups, if these are named.

- compartments:

  A vector for the model compartment names.

## Value

A `<data.table>` with the columns "compartment", "demography_group",
"value", and "time"; these specify the epidemiological compartment, the
name of the demography group, the number of individuals of that group in
the compartment, and the model timestep, respectively. Names for the
demographic groups are generated if no names are provided in the
`population` object; these are of the form "demo_group_X".
