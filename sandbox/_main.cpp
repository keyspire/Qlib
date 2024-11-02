#include "_main.h"

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

auto table_student_t_95p = std::vector<double>{
    12.7062047362, 4.30265272975, 3.18244630528, 2.77644510520, 2.57058183564,
    2.44691185114, 2.36462425159, 2.30600413520, 2.26215716280, 2.22813885199,
    2.20098516009, 2.17881282967, 2.16036865646, 2.14478668792, 2.13144954556,
    2.11990529922, 2.10981557783, 2.10092204024, 2.09302405441, 2.08596344727,
    2.07961384473, 2.07387306790, 2.06865761042, 2.06389856163, 2.05953855275,
    2.05552943864, 2.05183051648, 2.04840714180, 2.04522964213, 2.04227245630};

analyze::feed_normal_rv::feed_normal_rv(uint64_t seed)
    : _engine{std::mt19937_64{seed}},
      _dist{std::normal_distribution<>{0., 1.}} {}

double analyze::feed_normal_rv::get_single_normal_rv() {
  return _dist(_engine);
}

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

std::vector<analyze::result_t> analyze::get_95p_confidence_interval_t_test(
    const std::vector<result_t> divided_results) {
  auto _size = divided_results.size();
  auto _variance_vec = analyze::result_t{};
  for (auto i = 0; i < _size; ++i) {
    if (divided_results.at(i).size() > 30 || divided_results.at(i).size() < 1) {
      throw std::logic_error(
          "we have not implemented t-distribution value yet.");
    }
    _variance_vec.push_back(analyze::variance(divided_results.at(i)));
  }

  auto result = std::vector<analyze::result_t>{};
  for (auto i = 0; i < _size; ++i) {
    auto _mean = std::reduce(divided_results.at(i).begin(),
                             divided_results.at(i).end(), 0.0);
    auto _n = divided_results.at(i).size();
    _mean /= (static_cast<double>(_n));
    auto _variance = _variance_vec.at(i);
    double _r_minus;
    double _r_plus;
    _r_minus =
        _mean - table_student_t_95p.at(_n - 1) * std::sqrt(_variance / _n);
    _r_plus =
        _mean + table_student_t_95p.at(_n - 1) * std::sqrt(_variance / _n);
    result.push_back(analyze::result_t{_r_minus, _r_plus});
  }

  return result;
}

// åpè≥ÇÕÇµÇ»Ç¢ÇŸÇ§Ç™Ç¢Ç¢ÇÁÇµÇ¢
int main() {
  auto random_0{analyze::feed_normal_rv{428}};
  auto random_5489{analyze::feed_normal_rv{5489}};

  constexpr uint64_t num_runs = 50000;

  double log_current_val_0;
  double log_current_val_5489;
  double log_current_val_0_neg;
  double log_current_val_5489_neg;
  double fx_rate;

  constexpr double spot = 100.0;
  const double log_spot = std::log(spot);
  constexpr double dt = 1.;
  constexpr double vol = 0.2;
  constexpr double rate_foreign = 0.;
  constexpr double rate_base = 0.;
  const double numeraire_rate = rate_base;
  constexpr double rfr = rate_base - rate_foreign;
  constexpr double strike = 100.0;

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

    vec_payoff_0.push_back(payoff_0);
    vec_payoff_5489.push_back(payoff_5489);
  }

  const double discount_f = std::exp(-numeraire_rate * 1.0);

  double analytic =
      european_fx_call_value(dt, spot, strike, vol, rate_foreign, rate_base);

  double analytic_fxrate_expect =
      expectation_fxrate_at_maturity(dt, spot, rate_foreign, rate_base);

  analyze::result_t vec_payoff_0_mod;
  analyze::result_t vec_payoff_5489_mod;
  uint64_t _sum = 1;
  for (auto i = 0; i < num_runs; ++i) {
    vec_payoff_0_mod.push_back(vec_payoff_0.at(i) * discount_f / (2. * _sum));
    vec_payoff_5489_mod.push_back(vec_payoff_5489.at(i) * discount_f /
                                  (2. * _sum));
    _sum += 1;
  }

  auto error_bar_0 = analyze::get_95p_confidence_interval_t_test(
      analyze::get_partition(vec_payoff_0_mod, 30));
  auto error_bar_5489 = analyze::get_95p_confidence_interval_t_test(
      analyze::get_partition(vec_payoff_5489_mod, 30));

  // console out
  std::cout << "base_rate: " << rate_base << ", foreign_rate: " << rate_foreign
            << ", spot = " << spot << ", strike = " << strike
            << ", vol: " << vol << ", term (y):" << dt << std::endl;
  std::cout << "mc: negative_corr method, num_run: " << num_runs << std::endl
            << "seed = 0 result: " << payoff_0 * discount_f / (2. * num_runs)
            << std::endl
            << "95% interval lower side: " << error_bar_0.back().at(0)
            << ", 95% interval higher side: " << error_bar_0.back().at(1)
            << std::endl

            << "seed = 5489 result: "
            << payoff_5489 * discount_f / (2. * num_runs) << std::endl
            << "95% interval lower side: " << error_bar_5489.back().at(0)
            << ", 95% interval higher side: " << error_bar_5489.back().at(1)
            << std::endl

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
