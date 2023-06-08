# Printing intervention class

    Code
      close_schools
    Output
      <intervention>
      Intervention name: "close_schools"
      Time begin:
      [1] 100
      Time end:
      [1] 150
      Contact reduction:
           [,1]
      [1,]  0.2
      [2,]  0.0

# Concatenating `intervention`s works

    Code
      multi_npi
    Output
      <intervention>
      Intervention name: NA
      Time begin:
      npi_1 npi_2 
         30    45 
      Time end:
      npi_1 npi_2 
         60    75 
      Contact reduction:
           npi_1 npi_2
      [1,]  0.15   0.1
      [2,]  0.15   0.1
      [3,]  0.15   0.1

