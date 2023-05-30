# Printing vaccination class

    Code
      elder_vaccination
    Output
      <vaccination>
      Vaccination name: "elder_vaccination"
      
      Time begin:
           dose_1
      [1,]      0
      [2,]      0
      [3,]      0
      
      Time end:
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
    Output
      <vaccination>
      Vaccination name: "triple_vaccination"
      
      Time begin:
           dose_1 dose_2 dose_3
      [1,]      0      0      0
      
      Time end:
           dose_1 dose_2 dose_3
      [1,]     31     31     31
      
      Vaccination rate:
           dose_1 dose_2 dose_3
      [1,]  1e-04  1e-04  1e-04

# Multi-dose vaccination using `c()`

    Code
      c(vax_1, vax_2)
    Output
      <vaccination>
      Vaccination name: "vax_regime"
      
      Time begin:
           dose_1 dose_2
      [1,]      1    101
      
      Time end:
           dose_1 dose_2
      [1,]    100    200
      
      Vaccination rate:
           dose_1 dose_2
      [1,]  0.001  0.001

