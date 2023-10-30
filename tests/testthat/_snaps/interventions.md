# Printing <contacts_intervention> class

    Code
      close_schools
    Message
      <contacts_intervention> object
    Output
      
       Intervention name: 
    Message
      "close_schools"
    Output
      
       Begins at: 
      [1] 100
      
       Ends at: 
      [1] 150
      
       Reduction: 
                   Interv. 1
      Demo. grp. 1       0.2
      Demo. grp. 2       0.0

# Concatenating `intervention`s works

    Code
      multi_npi
    Message
      <contacts_intervention> object
    Output
      
       Intervention name: 
    Message
      NA
    Output
      
       Begins at: 
           npi_1 npi_2
      [1,]    30    45
      
       Ends at: 
           npi_1 npi_2
      [1,]    60    75
      
       Reduction: 
                   Interv. 1 Interv. 2
      Demo. grp. 1      0.15       0.1
      Demo. grp. 2      0.15       0.1
      Demo. grp. 3      0.15       0.1

