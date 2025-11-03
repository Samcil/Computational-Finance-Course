#include <Rcpp.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;
using namespace RcppParallel;

namespace {

class DeltaWorker : public Worker {
public:
  DeltaWorker(const NumericVector& terminal_prices,
              const NumericVector& strikes,
              double initial_price,
              double discount_factor,
              double sign_multiplier,
              bool is_call,
              NumericVector& output)
    : terminal_prices_(terminal_prices),
      strikes_(strikes),
      initial_price_(initial_price),
      discount_factor_(discount_factor),
      sign_multiplier_(sign_multiplier),
      is_call_(is_call),
      output_(output) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n_prices = static_cast<std::size_t>(terminal_prices_.size());
    if (n_prices == 0) {
      stop("`terminal_prices` must contain at least one observation.");
    }

    const double scale = discount_factor_ * sign_multiplier_ /
      (initial_price_ * static_cast<double>(n_prices));

    for (std::size_t idx = begin; idx < end; ++idx) {
      const double strike = strikes_[idx];
      double sum = 0.0;

      for (std::size_t price_idx = 0; price_idx < n_prices; ++price_idx) {
        const double price = terminal_prices_[price_idx];
        const bool indicator = is_call_ ? (price > strike) : (price < strike);
        if (indicator) {
          sum += price;
        }
      }

      output_[idx] = scale * sum;
    }
  }

private:
  const RVector<double> terminal_prices_;
  const RVector<double> strikes_;
  const double initial_price_;
  const double discount_factor_;
  const double sign_multiplier_;
  const bool is_call_;
  RVector<double> output_;
};

class VegaWorker : public Worker {
public:
  VegaWorker(const NumericVector& terminal_prices,
             const NumericVector& adjustment,
             const NumericVector& strikes,
             double volatility,
             double discount_factor,
             double sign_multiplier,
             bool is_call,
             NumericVector& output)
    : terminal_prices_(terminal_prices),
      adjustment_(adjustment),
      strikes_(strikes),
      volatility_(volatility),
      discount_factor_(discount_factor),
      sign_multiplier_(sign_multiplier),
      is_call_(is_call),
      output_(output) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n_prices = static_cast<std::size_t>(terminal_prices_.size());
    if (n_prices == 0) {
      stop("`terminal_prices` must contain at least one observation.");
    }
    if (adjustment_.size() != static_cast<int>(n_prices)) {
      stop("`adjustment` must have the same length as `terminal_prices`.");
    }

    const double scale = discount_factor_ * sign_multiplier_ /
      (volatility_ * static_cast<double>(n_prices));

    for (std::size_t idx = begin; idx < end; ++idx) {
      const double strike = strikes_[idx];
      double sum = 0.0;

      for (std::size_t price_idx = 0; price_idx < n_prices; ++price_idx) {
        const double price = terminal_prices_[price_idx];
        const bool indicator = is_call_ ? (price > strike) : (price < strike);
        if (indicator) {
          sum += price * adjustment_[price_idx];
        }
      }

      output_[idx] = scale * sum;
    }
  }

private:
  const RVector<double> terminal_prices_;
  const RVector<double> adjustment_;
  const RVector<double> strikes_;
  const double volatility_;
  const double discount_factor_;
  const double sign_multiplier_;
  const bool is_call_;
  RVector<double> output_;
};

} // namespace

// [[Rcpp::export]]
NumericVector pathwise_delta_cpp(NumericVector terminal_prices,
                                 NumericVector strikes,
                                 double initial_price,
                                 double discount_factor,
                                 double sign_multiplier,
                                 bool is_call) {
  if (initial_price <= 0) {
    stop("`initial_price` must be positive.");
  }
  if (discount_factor < 0) {
    stop("`discount_factor` must be non-negative.");
  }
  if (strikes.size() == 0) {
    stop("`strikes` must contain at least one element.");
  }

  NumericVector result(strikes.size());
  DeltaWorker worker(
    terminal_prices,
    strikes,
    initial_price,
    discount_factor,
    sign_multiplier,
    is_call,
    result
  );
  parallelFor(0, result.size(), worker);
  return result;
}

// [[Rcpp::export]]
NumericVector pathwise_vega_cpp(NumericVector terminal_prices,
                                NumericVector adjustment,
                                NumericVector strikes,
                                double volatility,
                                double discount_factor,
                                double sign_multiplier,
                                bool is_call) {
  if (volatility <= 0) {
    stop("`volatility` must be positive.");
  }
  if (discount_factor < 0) {
    stop("`discount_factor` must be non-negative.");
  }
  if (strikes.size() == 0) {
    stop("`strikes` must contain at least one element.");
  }

  NumericVector result(strikes.size());
  VegaWorker worker(
    terminal_prices,
    adjustment,
    strikes,
    volatility,
    discount_factor,
    sign_multiplier,
    is_call,
    result
  );
  parallelFor(0, result.size(), worker);
  return result;
}
