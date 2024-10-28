#pragma once
#include <memory>
#include <vector>

#include "../qpricer_lib_core/montecarlo.hpp"
#include "./template/template.h"

namespace qlib::bs::mc {
namespace mc = qlib::core::mc;

struct bs_update_setting final : mc::random_process_update_base_setting {
  bs_update_setting() = delete;
  bs_update_setting(double rfr, const std::vector<double>& gearing_1,
                    const std::vector<double>& gearing_2,
                    const std::vector<double>& volatilities,
                    const double initial_val,
                    const std::shared_ptr<mc::random_system_setting> rsconfig,
                    const uint_fast64_t index_in_system);
  const std::vector<double> gearing_1;
  const std::vector<double> gearing_2;
  const double rfr;
  const std::vector<double> volatilities;
  const uint_fast64_t steps;
};

class bs_random_process final : mc::i_random_process {
 public:
  bs_random_process() = delete;
  bs_random_process(const bs_update_setting& setting);

  mc::path_t get_path(
      const std::vector<std::vector<double>>& uncertainties_increments,
      const std::vector<double>& time_increments,
      const uint_fast64_t steps) override;

 private:
  void _update(const std::vector<double>& uncertainties_increment,
               const double time_increment,
               const uint_fast64_t current_step) override;
  bs_update_setting _setting;
  double current_val;
  mc::path_t path;
};

template <class... BsRps>
class bs_random_system : mc::i_random_system<BsRps...> {
  static_assert(qlib::utils::core::all_base_of<bs_random_process, BsRps...>(),
                "bs_random_system is a class whose random processes are only "
                "of Black-Scholes type.");

 public:
  bs_random_system() = delete;
  bs_random_system(const std::tuple<BsRps...>& rps,
                   const mc::random_system_setting& setting);
  std::unique_ptr<std::vector<mc::path_t>> get_paths() override;
  uint_fast64_t get_last_seed();

 private:
  template <int N>
  void _register_path_for_each_tuples();
  std::tuple<BsRps...> _rps;
  const mc::random_system_setting _setting;
  uint_fast64_t _seed;
  // mc::relation_vectors_t _rel_mat_trimmed;
  std::vector<std::vector<double>> _uncertainties_increments;
  std::vector<mc::path_t> _paths;
  std::normal_distribution<> _dist{0.0, 1.0};
};

}  // namespace qlib::bs::mc