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
      [0,40): 31,325,592 (50%)
      [40,65): 19,288,101 (30%)
      65+: 9,673,058 (20%)
      
       Contact matrix 
                       
      contact.age.group    [0,40)  [40,65)      65+
                [0,40)  9.6963190 4.792619 2.235873
                [40,65) 2.9509586 4.178218 2.468756
                65+     0.6904172 1.238091 1.714286

---

    structure(list(name = "UK population", contact_matrix = structure(c(9.69631901840491, 
    2.95095856108701, 0.690417159802703, 4.79261923677809, 4.17821782178218, 
    1.23809069406646, 2.23587269483738, 2.46875583236594, 1.71428571428571
    ), .Dim = c(3L, 3L), .Dimnames = list(contact.age.group = NULL, 
        NULL)), demography_vector = c(31325592, 19288101, 9673058
    ), initial_conditions = structure(c(0.9999, 0.9999, 0.9999, 5e-05, 
    5e-05, 5e-05, 5e-05, 5e-05, 5e-05, 0, 0, 0, 0, 0, 0), .Dim = c(3L, 
    5L))), class = "population")

