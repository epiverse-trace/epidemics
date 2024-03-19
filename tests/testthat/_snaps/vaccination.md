# Printing vaccination class

    Code
      elder_vaccination
    Message
      <vaccination> object
    Output
      
       Vaccination name: 
    Message
      "elder_vaccination"
    Output
      
       Begins at: 
           dose_1
      [1,]      0
      [2,]      0
      [3,]      0
      
       Ends at: 
           dose_1
      [1,]    200
      [2,]    200
      [3,]    200
      
       Vaccination rate: 
           dose_1
      [1,]  0e+00
      [2,]  0e+00
      [3,]  1e-04

---

    Code
      triple_vaccination
    Message
      <vaccination> object
    Output
      
       Vaccination name: 
    Message
      "triple_vaccination"
    Output
      
       Begins at: 
           dose_1 dose_2 dose_3
      [1,]      0      0      0
      
       Ends at: 
           dose_1 dose_2 dose_3
      [1,]     31     31     31
      
       Vaccination rate: 
           dose_1 dose_2 dose_3
      [1,]  1e-04  1e-04  1e-04

# Multi-dose vaccination using `c()`

    Code
      c(vax_1, vax_2)
    Message
      <vaccination> object
    Output
      
       Vaccination name: 
    Message
      "vax_regime"
    Output
      
       Begins at: 
           dose_1 dose_2
      [1,]      1    101
      
       Ends at: 
           dose_1 dose_2
      [1,]    100    200
      
       Vaccination rate: 
           dose_1 dose_2
      [1,]  0.001  0.001

# Epidemic model with vaccination

    Code
      data_vaccination[data_vaccination$compartment == "vaccinated" &
        data_vaccination$time < 5, ]
    Output
           time demography_group compartment     value
          <num>           <char>      <char>     <num>
       1:     0           [0,40)  vaccinated    0.0000
       2:     0          [40,65)  vaccinated    0.0000
       3:     0              65+  vaccinated    0.0000
       4:     1           [0,40)  vaccinated    0.0000
       5:     1          [40,65)  vaccinated    0.0000
       6:     1              65+  vaccinated  967.3058
       7:     2           [0,40)  vaccinated    0.0000
       8:     2          [40,65)  vaccinated    0.0000
       9:     2              65+  vaccinated 1934.6116
      10:     3           [0,40)  vaccinated    0.0000
      11:     3          [40,65)  vaccinated    0.0000
      12:     3              65+  vaccinated 2901.9174
      13:     4           [0,40)  vaccinated    0.0000
      14:     4          [40,65)  vaccinated    0.0000
      15:     4              65+  vaccinated 3869.2232

