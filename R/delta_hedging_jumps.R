#' Simulate Delta Hedging Strategy with Jump Risk
#'
#' Evaluates Black-Scholes style delta hedging when the underlying follows a
#' jump-diffusion process (e.g. Merton). The hedging logic mirrors
#' [simulate_delta_hedge_bs()] but accepts paths generated from jump models so
#' that residual jump risk is captured in the hedging error.
#'
#' @param spot_paths Tibble of simulated paths with columns `path_id`, `time`,
#'   and `stock_price`.
#' @param option_spec Black-Scholes option specification describing the payoff.
#' @param diffusion_volatility Numeric. Diffusive volatility assumed by the
#'   hedging strategy. This is typically the instantaneous volatility parameter
#'   of the jump-diffusion model.
#' @param hedge_steps Integer. Number of trading dates between 0 and maturity.
#' @param transaction_cost Numeric proportional cost per trade notional.
#' @param initial_capital Optional starting capital. Defaults to the
#'   Black-Scholes price implied by `diffusion_volatility`.
#'
#' @return A tibble identical in structure to the output of
#'   [simulate_delta_hedge_bs()].
#'
#' @export
simulate_delta_hedge_jumps <- function(spot_paths,
																			 option_spec,
																			 diffusion_volatility,
																			 hedge_steps,
																			 transaction_cost = 0,
																			 initial_capital = NULL) {
	simulate_delta_hedge_bs(
		spot_paths = spot_paths,
		option_spec = option_spec,
		volatility = diffusion_volatility,
		hedge_steps = hedge_steps,
		transaction_cost = transaction_cost,
		initial_capital = initial_capital
	)
}
