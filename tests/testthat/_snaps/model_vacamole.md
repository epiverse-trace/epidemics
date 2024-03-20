# Vacamole model: two dose vaccination and stats. correctness

    Code
      data[grepl("dose", data$compartment, fixed = TRUE) & time %in% seq(50, 55), ]
    Output
           time demography_group         compartment       value
          <num>           <char>              <char>       <num>
       1:    50           [0,60) vaccinated_one_dose 2370379.149
       2:    50              60+ vaccinated_one_dose  638905.568
       3:    50           [0,60) vaccinated_two_dose    3957.270
       4:    50              60+ vaccinated_two_dose    1066.626
       5:    51           [0,60) vaccinated_one_dose 2394121.139
       6:    51              60+ vaccinated_one_dose  645305.113
       7:    51           [0,60) vaccinated_two_dose   27700.882
       8:    51              60+ vaccinated_two_dose    7466.380
       9:    52           [0,60) vaccinated_one_dose 2417863.063
      10:    52              60+ vaccinated_one_dose  651704.650
      11:    52           [0,60) vaccinated_two_dose   51444.480
      12:    52              60+ vaccinated_two_dose   13866.132
      13:    53           [0,60) vaccinated_one_dose 2441604.917
      14:    53              60+ vaccinated_one_dose  658104.178
      15:    53           [0,60) vaccinated_two_dose   75188.064
      16:    53              60+ vaccinated_two_dose   20265.882
      17:    54           [0,60) vaccinated_one_dose 2465346.699
      18:    54              60+ vaccinated_one_dose  664503.696
      19:    54           [0,60) vaccinated_two_dose   98931.632
      20:    54              60+ vaccinated_two_dose   26665.630
      21:    55           [0,60) vaccinated_one_dose 2489088.407
      22:    55              60+ vaccinated_one_dose  670903.205
      23:    55           [0,60) vaccinated_two_dose  122675.184
      24:    55              60+ vaccinated_two_dose   33065.377
           time demography_group         compartment       value

