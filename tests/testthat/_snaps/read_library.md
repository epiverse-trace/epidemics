# Reading model library JSON

    Code
      ml
    Output
        model_type model_name         model_function            model_args_checker
      1   epidemic    default  .epidemic_default_cpp  .check_args_epidemic_default
      2   epidemic      ebola    .epidemic_ebola_cpp    .check_args_epidemic_ebola
      3   epidemic   vacamole .epidemic_vacamole_cpp .check_args_epidemic_vacamole
                     model_args_prepper
      1  .prepare_args_epidemic_default
      2    .prepare_args_epidemic_ebola
      3 .prepare_args_epidemic_vacamole
                                                                                                                                                                         compartments
      1                                                                                                                       susceptible, exposed, infectious, recovered, vaccinated
      2                                                                                                                                   susceptible, exposed, infectious, recovered
      3 susceptible, vaccinated_one_dose, vaccinated_two_dose, exposed, exposed_vaccinated, infectious, infectious_vaccinated, hospitalised, hospitalised_vaccinated, dead, recovered

