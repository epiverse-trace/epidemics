# Vacamole model: two dose vaccination and stats. correctness

    Code
      data[grepl("dose", data$compartment, fixed = TRUE) & time %in% seq(50, 55), ]
    Output
           time demography_group         compartment        value
          <num>           <char>              <char>        <num>
       1:    50           [0,60) vaccinated_one_dose 2.408161e+06
       2:    50         [60,Inf) vaccinated_one_dose 6.559418e+05
       5:    51           [0,60) vaccinated_one_dose 2.432242e+06
       6:    51         [60,Inf) vaccinated_one_dose 6.625010e+05
       7:    51           [0,60) vaccinated_two_dose 2.408187e+04
       8:    51         [60,Inf) vaccinated_two_dose 6.559451e+03
       9:    52           [0,60) vaccinated_one_dose 2.456322e+06
      10:    52         [60,Inf) vaccinated_one_dose 6.690602e+05
      11:    52           [0,60) vaccinated_two_dose 4.816372e+04
      12:    52         [60,Inf) vaccinated_two_dose 1.311890e+04
      13:    53           [0,60) vaccinated_one_dose 2.480402e+06
      14:    53         [60,Inf) vaccinated_one_dose 6.756195e+05
      15:    53           [0,60) vaccinated_two_dose 7.224556e+04
      16:    53         [60,Inf) vaccinated_two_dose 1.967835e+04
      17:    54           [0,60) vaccinated_one_dose 2.504482e+06
      18:    54         [60,Inf) vaccinated_one_dose 6.821787e+05
      19:    54           [0,60) vaccinated_two_dose 9.632738e+04
      20:    54         [60,Inf) vaccinated_two_dose 2.623779e+04
      21:    55           [0,60) vaccinated_one_dose 2.528562e+06
      22:    55         [60,Inf) vaccinated_one_dose 6.887379e+05
      23:    55           [0,60) vaccinated_two_dose 1.204092e+05
      24:    55         [60,Inf) vaccinated_two_dose 3.279723e+04
           time demography_group         compartment        value
          <num>           <char>              <char>        <num>

