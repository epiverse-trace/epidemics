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
                       age.group
      contact.age.group    [0,40)  [40,65) [65,Inf)
               [0,40)   9.6963190 4.784327 2.240342
               [40,65)  2.9570904 4.178218 2.478450
               [65,Inf) 0.6863798 1.228532 1.714286
      
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
      Dem.grp.1:31,844,217(50%)
      Dem.grp.2:19,682,231(30%)
      Dem.grp.3:9,756,203(20%)
      
      Contactmatrix

      age.group
      contact.age.group[0,40)[40,65)[65,Inf)
      Dem.grp.1:9.69631904.7843272.240342
      Dem.grp.2:2.95709044.1782182.478450
      Dem.grp.3:0.68637981.2285321.714286
      
       Initial Conditions 
             [,1]  [,2]  [,3] [,4] [,5]
      [1,] 0.9999 5e-05 5e-05    0    0
      [2,] 0.9999 5e-05 5e-05    0    0
      [3,] 0.9999 5e-05 5e-05    0    0

