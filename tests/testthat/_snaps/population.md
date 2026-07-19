# Printing population class

    Code
      uk_population
    Message
      <population> object
    Output
      
       Population name: 
    Message
      "UK population"
    Output
      
       Demography 
      [0,40): 31,844,217 (50%)
      [40,65): 19,682,231 (30%)
      [65,Inf): 9,756,203 (20%)
      
       Contact matrix 
                contact.age.group
      age.group    [0,40)  [40,65)  [65,Inf)
        [0,40)   9.696319 2.957090 0.6863798
        [40,65)  4.784327 4.178218 1.2285324
        [65,Inf) 2.240342 2.478450 1.7142857
      
       Initial Conditions 
             [,1]  [,2]  [,3] [,4] [,5]
      [1,] 0.9999 5e-05 5e-05    0    0
      [2,] 0.9999 5e-05 5e-05    0    0
      [3,] 0.9999 5e-05 5e-05    0    0

# Printing population class with no row or column labels

    Code
      uk_population
    Message
      <population> object
    Output
      
       Population name: 
    Message
      "UK population"
    Output
      
       Demography 
      Dem. grp. 1: 31,844,217 (50%)
      Dem. grp. 2: 19,682,231 (30%)
      Dem. grp. 3: 9,756,203 (20%)
      
       Contact matrix 
                    contact.age.group
      age.group        [0,40)  [40,65)  [65,Inf)
        Dem. grp. 1: 9.696319 2.957090 0.6863798
        Dem. grp. 2: 4.784327 4.178218 1.2285324
        Dem. grp. 3: 2.240342 2.478450 1.7142857
      
       Initial Conditions 
             [,1]  [,2]  [,3] [,4] [,5]
      [1,] 0.9999 5e-05 5e-05    0    0
      [2,] 0.9999 5e-05 5e-05    0    0
      [3,] 0.9999 5e-05 5e-05    0    0

