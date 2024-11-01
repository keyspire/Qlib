#include "_main.h"

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

feed_normal_rv::feed_normal_rv(uint64_t seed)
    : _engine{std::mt19937_64{seed}},
      _dist{std::normal_distribution<>{0., 1.}} {}

double feed_normal_rv::get_single_normal_rv() { return _dist(_engine); }

double update_bs_process(double log_current_val, const double dt,
                         const double vol, const double rfr,
                         const double standard_normal_rv) {
  log_current_val +=
      (rfr - 0.5 * vol * vol) * dt + vol * std::sqrt(dt) * standard_normal_rv;
  return log_current_val;
}

// vol means stdev..
double d_plus(const double term, const double spot, const double strike,
              const double vol, const double rfr) {
  return (std::log(spot / strike) + (rfr + 0.5 * vol * vol) * term) /
         (vol * std::sqrt(term));
}
double d_minus(const double term, const double spot, const double strike,
               const double vol, const double rfr) {
  return (std::log(spot / strike) + (rfr - 0.5 * vol * vol) * term) /
         (vol * std::sqrt(term));
}

double cumulative_normal_dist(const double x, const double mu,
                              const double vol) {
  return 0.5 * (1 + std::erf((x - mu) / (std::sqrt(2) * vol)));
}

double european_fx_call_value(const double term, const double spot,
                              const double strike, const double vol,
                              const double rate_foreign,
                              const double rate_base) {
  return std::exp(-rate_foreign * term) * spot *
             cumulative_normal_dist(
                 d_plus(term, spot, strike, vol, rate_base - rate_foreign), 0,
                 1) -
         std::exp(-rate_base * term) * strike *
             cumulative_normal_dist(
                 d_minus(term, spot, strike, vol, rate_base - rate_foreign), 0,
                 1);
}

double expectation_fxrate_at_maturity(const double term, const double spot,
                                      const double rate_foreign,
                                      const double rate_base) {
  return spot * std::exp((rate_base - rate_foreign) * term);
}

// åpè≥ÇÕÇµÇ»Ç¢ÇŸÇ§Ç™Ç¢Ç¢ÇÁÇµÇ¢
int main() {
  auto random_0{feed_normal_rv{0}};
  auto random_5489{feed_normal_rv{5489}};

  constexpr double spot = 100;
  constexpr uint64_t num_runs = 500000;

  const double log_spot = std::log(spot);

  double log_current_val_0;
  double log_current_val_5489;
  double log_current_val_0_neg;
  double log_current_val_5489_neg;
  double fx_rate;

  constexpr double dt = 1.;
  constexpr double vol = 0.25;
  constexpr double numeraire_rate = 0.01;
  constexpr double rate_foreign = 0.05;
  constexpr double rate_base = 0.01;
  constexpr double rfr = -0.04;
  constexpr double strike = 100;

  double payoff_0{0.};
  double payoff_5489{0.};

  double payoff_0_cv{0.};
  double payoff_5489_cv{0.};

  double coeff{0.5};

  std::vector<double> vec_payoff_0;
  std::vector<double> vec_payoff_5489;

  for (uint64_t i = 0; i < num_runs; ++i) {
    auto realized_0 = random_0.get_single_normal_rv();
    auto realized_5489 = random_5489.get_single_normal_rv();

    log_current_val_0 = log_spot;
    log_current_val_5489 = log_spot;
    log_current_val_0_neg = log_spot;
    log_current_val_5489_neg = log_spot;

    fx_rate = spot;

    log_current_val_0 =
        update_bs_process(log_current_val_0, dt, vol, rfr, realized_0);
    log_current_val_5489 =
        update_bs_process(log_current_val_5489, dt, vol, rfr, realized_5489);

    log_current_val_0_neg =
        update_bs_process(log_current_val_0_neg, dt, vol, rfr, -realized_0);
    log_current_val_5489_neg = update_bs_process(log_current_val_5489_neg, dt,
                                                 vol, rfr, -realized_5489);

    payoff_0 += std::max(std::exp(log_current_val_0) - strike, 0.) +
                std::max(std::exp(log_current_val_0_neg) - strike, 0.);
    payoff_5489 += std::max(std::exp(log_current_val_5489) - strike, 0.) +
                   std::max(std::exp(log_current_val_5489_neg) - strike, 0.);

    payoff_0_cv += std::max(std::exp(log_current_val_0) - strike, 0.) +
                   std::max(std::exp(log_current_val_0_neg) - strike, 0.) -
                   coeff * std::exp(log_current_val_0) -
                   coeff * std::exp(log_current_val_0_neg);
    payoff_5489_cv +=
        std::max(std::exp(log_current_val_5489) - strike, 0.) +
        std::max(std::exp(log_current_val_5489_neg) - strike, 0.) -
        coeff * std::exp(log_current_val_5489) -
        coeff * std::exp(log_current_val_5489_neg);

    /* vec_payoff_0.push_back(payoff_0);
     vec_payoff_5489.push_back(payoff_5489);*/
  }

  const double discount_f = std::exp(-numeraire_rate * 1.0);

  double analytic =
      european_fx_call_value(dt, spot, strike, vol, rate_foreign, rate_base);

  double analytic_fxrate_expect =
      expectation_fxrate_at_maturity(dt, spot, rate_foreign, rate_base);

  std::cout << "base_rate: " << rate_base << ", foreign_rate: " << rate_foreign
            << ", ATM spot = 100, vol: " << vol << ", term (y):" << dt
            << std::endl;
  std::cout << "mc: negative_corr method, num_run: " << num_runs << std::endl
            << "seed = 0 result: " << payoff_0 * discount_f / (2. * num_runs)
            << std::endl

            << "seed = 5489 result: "
            << payoff_5489 * discount_f / (2. * num_runs) << std::endl
            << "mc: negative_corr + control_val method, num_run: " << num_runs
            << std::endl
            << "seed = 0 result: "
            << discount_f * (payoff_0_cv / (2. * num_runs) +
                             coeff * analytic_fxrate_expect)

            << std::endl

            << "seed = 5489 result: "
            << discount_f * (payoff_5489_cv / (2. * num_runs) +
                             coeff * analytic_fxrate_expect)
            << std::endl

            << "analytic: " << analytic << std::endl;
  return 0;
}