# Printing <contacts_intervention> class

    Code
      close_schools
    Output
      <contacts_intervention>
      Intervention name: "close_schools"
      
      Time begin:
      [1] 100
      
      Time end:
      [1] 150
      
      Reduction:
                   Interv. 1
      Demo. grp. 1       0.2
      Demo. grp. 2       0.0

# Concatenating `intervention`s works

    Code
      multi_npi
    Output
      <contacts_intervention>
      Intervention name: NA
      
      Time begin:
           npi_1 npi_2
      [1,]    30    45
      
      Time end:
           npi_1 npi_2
      [1,]    60    75
      
      Reduction:
                   Interv. 1 Interv. 2
      Demo. grp. 1      0.15       0.1
      Demo. grp. 2      0.15       0.1
      Demo. grp. 3      0.15       0.1

