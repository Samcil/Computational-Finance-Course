#include <Rcpp.h>
#include <RcppParallel.h>
#include <complex>

// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;

namespace {

class ThetaIntegrator : public Worker {
public:
  ThetaIntegrator(const NumericVector& theta,
                  const NumericVector& exp_neg,
                  const NumericVector& u_vals,
                  double step_size,
                  double lambda_val,
                  NumericVector& real_out,
                  NumericVector& imag_out)
    : theta_(theta),
      exp_neg_(exp_neg),
      u_vals_(u_vals),
      step_size_(step_size),
      lambda_(lambda_val),
      real_out_(real_out),
      imag_out_(imag_out) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n = theta_.size();

    for (std::size_t idx = begin; idx < end; ++idx) {
      const double u = u_vals_[idx];
      const std::complex<double> factor(-1.0 / lambda_, u / lambda_);
      std::complex<double> acc(0.0, 0.0);

      for (std::size_t i = 0; i < n; ++i) {
        const double trap_weight = (i == 0 || i == n - 1) ? 0.5 : 1.0;
        const double one_minus_exp = 1.0 - exp_neg_[i];
        const double scaled_theta = theta_[i] * one_minus_exp;
        acc += factor * scaled_theta * trap_weight;
      }

      acc *= step_size_;
      real_out_[idx] = acc.real();
      imag_out_[idx] = acc.imag();
    }
  }

private:
  const RVector<double> theta_;
  const RVector<double> exp_neg_;
  const RVector<double> u_vals_;
  const double step_size_;
  const double lambda_;
  RVector<double> real_out_;
  RVector<double> imag_out_;
};

} // namespace

// [[Rcpp::export]]
ComplexVector bshw_theta_integrals_cpp(NumericVector u,
                                       double lambda,
                                       NumericVector theta,
                                       NumericVector exp_neg,
                                       double step_size) {
  if (theta.size() != exp_neg.size()) {
    stop("`theta` and `exp_neg` must have matching lengths.");
  }
  if (theta.size() < 2) {
    stop("Integration grid must contain at least two points.");
  }
  if (lambda == 0.0) {
    stop("`lambda` must be non-zero.");
  }

  NumericVector real_part(u.size());
  NumericVector imag_part(u.size());

  ThetaIntegrator worker(theta, exp_neg, u, step_size, lambda, real_part, imag_part);
  parallelFor(0, u.size(), worker);

  ComplexVector result(u.size());
  for (R_xlen_t i = 0; i < u.size(); ++i) {
    Rcomplex value;
    value.r = real_part[i];
    value.i = imag_part[i];
    result[i] = value;
  }

  return result;
}
