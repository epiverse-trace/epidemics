# Vacamole model: two dose vaccination and stats. correctness

    Code
      data[grepl("dose", data$compartment, fixed = TRUE) & time %in% seq(50, 55), ]
    Output
           time demography_group         compartment        value
          <num>           <char>              <char>        <num>
       1:    50           [0,60) vaccinated_one_dose 2.374336e+06
       2:    50              60+ vaccinated_one_dose 6.399722e+05
       5:    51           [0,60) vaccinated_one_dose 2.398078e+06
       6:    51              60+ vaccinated_one_dose 6.463717e+05
       7:    51           [0,60) vaccinated_two_dose 2.374361e+04
       8:    51              60+ vaccinated_two_dose 6.399754e+03
       9:    52           [0,60) vaccinated_one_dose 2.421820e+06
      10:    52              60+ vaccinated_one_dose 6.527713e+05
      11:    52           [0,60) vaccinated_two_dose 4.748721e+04
      12:    52              60+ vaccinated_two_dose 1.279951e+04
      13:    53           [0,60) vaccinated_one_dose 2.445562e+06
      14:    53              60+ vaccinated_one_dose 6.591708e+05
      15:    53           [0,60) vaccinated_two_dose 7.123080e+04
      16:    53              60+ vaccinated_two_dose 1.919926e+04
      17:    54           [0,60) vaccinated_one_dose 2.469304e+06
      18:    54              60+ vaccinated_one_dose 6.655703e+05
      19:    54           [0,60) vaccinated_two_dose 9.497437e+04
      20:    54              60+ vaccinated_two_dose 2.559901e+04
      21:    55           [0,60) vaccinated_one_dose 2.493046e+06
      22:    55              60+ vaccinated_one_dose 6.719698e+05
      23:    55           [0,60) vaccinated_two_dose 1.187179e+05
      24:    55              60+ vaccinated_two_dose 3.199875e+04
           time demography_group         compartment        value

