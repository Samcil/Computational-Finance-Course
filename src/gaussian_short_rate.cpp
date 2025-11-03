#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>

// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;

namespace {

class GaussianShortRateWorker : public Worker {
public:
  GaussianShortRateWorker(const NumericMatrix& shocks,
                          const NumericVector& theta,
                          const NumericVector& initial_rates,
                          double mean_reversion,
                          double dt,
                          NumericMatrix& rates,
                          NumericMatrix& discounts)
    : shocks_(shocks),
      theta_(theta),
      initial_rates_(initial_rates),
      mean_reversion_(mean_reversion),
      dt_(dt),
      rates_(rates),
      discounts_(discounts) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n_steps = static_cast<std::size_t>(theta_.size());

    for (std::size_t path_idx = begin; path_idx < end; ++path_idx) {
      double rate_prev = initial_rates_[path_idx];
      double discount_prev = 1.0;
      rates_(path_idx, 0) = rate_prev;
      discounts_(path_idx, 0) = discount_prev;

      for (std::size_t step = 0; step < n_steps; ++step) {
        const double drift = (theta_[step] - mean_reversion_ * rate_prev) * dt_;
        const double shock = shocks_(path_idx, step);
        const double rate_next = rate_prev + drift + shock;
        const double discount_next = discount_prev * std::exp(-0.5 * (rate_prev + rate_next) * dt_);

        rates_(path_idx, step + 1) = rate_next;
        discounts_(path_idx, step + 1) = discount_next;

        rate_prev = rate_next;
        discount_prev = discount_next;
      }
    }
  }

private:
  const RMatrix<double> shocks_;
  const RVector<double> theta_;
  const RVector<double> initial_rates_;
  const double mean_reversion_;
  const double dt_;
  RMatrix<double> rates_;
  RMatrix<double> discounts_;
};

} // namespace

// [[Rcpp::export]]
List gaussian_short_rate_paths_cpp(NumericMatrix shocks,
                                   NumericVector theta,
                                   NumericVector initial_rate,
                                   double mean_reversion,
                                   double dt) {
  const std::size_t n_paths = static_cast<std::size_t>(shocks.nrow());
  const std::size_t n_steps = static_cast<std::size_t>(shocks.ncol());

  if (theta.size() != static_cast<int>(n_steps)) {
    stop("`theta` length must match number of columns in `shocks`.");
  }

  if (!(initial_rate.size() == 1 || initial_rate.size() == static_cast<int>(n_paths))) {
    stop("`initial_rate` must have length 1 or match number of paths.");
  }

  NumericVector initial_rates(n_paths);
  if (initial_rate.size() == 1) {
    std::fill(initial_rates.begin(), initial_rates.end(), initial_rate[0]);
  } else {
    std::copy(initial_rate.begin(), initial_rate.end(), initial_rates.begin());
  }

  NumericMatrix rates(n_paths, n_steps + 1);
  NumericMatrix discounts(n_paths, n_steps + 1);

  GaussianShortRateWorker worker(shocks, theta, initial_rates, mean_reversion, dt, rates, discounts);
  parallelFor(0, n_paths, worker);

  return List::create(
    Named("rates") = rates,
    Named("discounts") = discounts
  );
}

namespace {

class GaussianMeasureSwitchWorker : public Worker {
public:
  GaussianMeasureSwitchWorker(const NumericMatrix& shocks,
                              const NumericVector& theta,
                              const NumericVector& initial_rates,
                              const NumericVector& lambda,
                              double mean_reversion,
                              double volatility,
                              double dt,
                              NumericMatrix& rates_p,
                              NumericMatrix& rates_q,
                              NumericMatrix& discounts_p,
                              NumericMatrix& discounts_q,
                              NumericMatrix& radon)
    : shocks_(shocks),
      theta_(theta),
      initial_rates_(initial_rates),
      lambda_(lambda),
      mean_reversion_(mean_reversion),
      volatility_(volatility),
      dt_(dt),
      sqrt_dt_(std::sqrt(dt)),
      rates_p_(rates_p),
      rates_q_(rates_q),
      discounts_p_(discounts_p),
      discounts_q_(discounts_q),
      radon_(radon) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n_steps = static_cast<std::size_t>(theta_.size());

    for (std::size_t path_idx = begin; path_idx < end; ++path_idx) {
      double rate_p = initial_rates_[path_idx];
      double rate_q = initial_rates_[path_idx];
      double discount_p = 1.0;
      double discount_q = 1.0;
      double log_rn = 0.0;

      rates_p_(path_idx, 0) = rate_p;
      rates_q_(path_idx, 0) = rate_q;
      discounts_p_(path_idx, 0) = discount_p;
      discounts_q_(path_idx, 0) = discount_q;
      radon_(path_idx, 0) = 1.0;

      for (std::size_t step = 0; step < n_steps; ++step) {
        const double theta_val = theta_[step];
        const double lambda_val = lambda_[step];
        const double shock_standard = shocks_(path_idx, step);
        const double dW_p = sqrt_dt_ * shock_standard;
        const double shock_term_p = volatility_ * dW_p;

        const double drift_p = (theta_val - mean_reversion_ * rate_p + volatility_ * lambda_val) * dt_;
        const double rate_next_p = rate_p + drift_p + shock_term_p;
        const double discount_next_p = discount_p * std::exp(-0.5 * (rate_p + rate_next_p) * dt_);

        const double dW_q = dW_p + lambda_val * dt_;
        const double shock_term_q = volatility_ * dW_q;
        const double drift_q = (theta_val - mean_reversion_ * rate_q) * dt_;
        const double rate_next_q = rate_q + drift_q + shock_term_q;
        const double discount_next_q = discount_q * std::exp(-0.5 * (rate_q + rate_next_q) * dt_);

        log_rn += -lambda_val * dW_p - 0.5 * lambda_val * lambda_val * dt_;

        rates_p_(path_idx, step + 1) = rate_next_p;
        rates_q_(path_idx, step + 1) = rate_next_q;
        discounts_p_(path_idx, step + 1) = discount_next_p;
        discounts_q_(path_idx, step + 1) = discount_next_q;
        radon_(path_idx, step + 1) = std::exp(log_rn);

        rate_p = rate_next_p;
        rate_q = rate_next_q;
        discount_p = discount_next_p;
        discount_q = discount_next_q;
      }
    }
  }

private:
  const RMatrix<double> shocks_;
  const RVector<double> theta_;
  const RVector<double> initial_rates_;
  const RVector<double> lambda_;
  const double mean_reversion_;
  const double volatility_;
  const double dt_;
  const double sqrt_dt_;
  RMatrix<double> rates_p_;
  RMatrix<double> rates_q_;
  RMatrix<double> discounts_p_;
  RMatrix<double> discounts_q_;
  RMatrix<double> radon_;
};

} // namespace

// [[Rcpp::export]]
List gaussian_measure_switch_cpp(NumericMatrix shocks,
                                 NumericVector theta,
                                 NumericVector initial_rate,
                                 double mean_reversion,
                                 double volatility,
                                 double dt,
                                 NumericVector lambda) {
  const std::size_t n_paths = static_cast<std::size_t>(shocks.nrow());
  const std::size_t n_steps = static_cast<std::size_t>(shocks.ncol());

  if (theta.size() != static_cast<int>(n_steps)) {
    stop("`theta` length must match number of columns in `shocks`.");
  }

  if (lambda.size() != static_cast<int>(n_steps)) {
    stop("`lambda` length must match number of columns in `shocks`.");
  }

  if (!(initial_rate.size() == 1 || initial_rate.size() == static_cast<int>(n_paths))) {
    stop("`initial_rate` must have length 1 or match number of paths.");
  }

  NumericVector initial_rates(n_paths);
  if (initial_rate.size() == 1) {
    std::fill(initial_rates.begin(), initial_rates.end(), initial_rate[0]);
  } else {
    std::copy(initial_rate.begin(), initial_rate.end(), initial_rates.begin());
  }

  NumericMatrix rates_p(n_paths, n_steps + 1);
  NumericMatrix rates_q(n_paths, n_steps + 1);
  NumericMatrix discounts_p(n_paths, n_steps + 1);
  NumericMatrix discounts_q(n_paths, n_steps + 1);
  NumericMatrix radon(n_paths, n_steps + 1);

  GaussianMeasureSwitchWorker worker(
    shocks,
    theta,
    initial_rates,
    lambda,
    mean_reversion,
    volatility,
    dt,
    rates_p,
    rates_q,
    discounts_p,
    discounts_q,
    radon
  );

  parallelFor(0, n_paths, worker);

  return List::create(
    Named("rates_p") = rates_p,
    Named("rates_q") = rates_q,
    Named("discounts_p") = discounts_p,
    Named("discounts_q") = discounts_q,
    Named("radon") = radon
  );
}
