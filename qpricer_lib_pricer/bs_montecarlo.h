#pragma once
#include <memory>
#include <vector>

#include "../qpricer_lib_core/montecarlo.hpp"

namespace qlib::bs::mc {
namespace mc = qlib::core::mc;

struct bs_update_setting final : mc::random_process_update_base_setting {
  bs_update_setting()=delete;
  bs_update_setting(int dim, double rfr, std::vector<double>& gearing_1,
                    std::vector<double>& gearing_2,
                    const std::vector<const std::vector<double>>& volatilities,
                    const double initial_val,
                    const std::shared_ptr<mc::random_system_setting> rsconfig);
  const std::vector<double> gearing_1;
  const std::vector<double> gearing_2;
  const double rfr;
  const std::vector<const std::vector<double>> volatilities;
};

class bs_random_process : mc::i_random_process<bs_update_setting> {
 public:
  bs_random_process() = delete;
  bs_random_process(const bs_update_setting& setting);

  mc::path_t get_path(
      const std::vector<const std::vector<double>>& uncertainties_increments,
      const std::vector<double>& time_increments) override;

 private:
  void _update(const std::vector<double> uncertainties_increment,
               const bs_update_setting& setting) override;
  bs_update_setting _setting;
  double current_val;
};
}  // namespace qlib::bs::mc