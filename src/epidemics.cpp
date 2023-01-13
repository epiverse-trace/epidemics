#include <Rcpp.h>
#include <sir_stochastic.h>

// [[Rcpp::interfaces(r, cpp)]]

//' @title Stochastic SIR epidemic simulation
//' @description Some description.
//' @param option A string giving the simulation option.
//' @param parameters A named list of parameters to the simulation.
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame epi_demic(const Rcpp::List &parameters, const Rcpp::String &option = "sir_stochastic") {
    if (option == "sir_stochastic") {
        return sir_stochastic(
            parameters["beta"], parameters["gamma"],
            parameters["N"], parameters["S0"], parameters["I0"],
            parameters["R0"], parameters["tf"]
        );
    } else {
        Rcpp::stop("Option not found, option must be 'sir_stochastic'");
    }
}
