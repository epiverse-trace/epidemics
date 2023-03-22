#include <Rcpp.h>
#include <epidemics.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]

//' @title An SIR model
//' @export
// [[Rcpp::export]]
Rcpp::List epidemic_default_cpp(const Rcpp::NumericVector &init,
                                const float &beta, const float &gamma) {
  // initial conditions from input
  odetools::state_type x = Rcpp::as<std::vector<double> >(init);

  // create a default epidemic with parameters
  epidemics::epidemic_default this_model(beta, gamma);

  //[ integrate_observ
  std::vector<odetools::state_type> x_vec;  // is a vector of double vectors
  std::vector<double> times;

  // a controlled stepper for constant step sizes
  boost::numeric::odeint::runge_kutta4<odetools::state_type> stepper;

  size_t steps = boost::numeric::odeint::integrate_const(
      stepper, this_model, x, 0.0, 200.0, 0.1,
      odetools::observer(x_vec, times));

  return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x_vec),
                            Rcpp::Named("time") = Rcpp::wrap(times));
  // return process_data(data);
}
