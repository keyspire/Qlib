#include "bs_montecarlo.h"

namespace qlib::bs::mc {
bs_update_setting::bs_update_setting(
    int dim, double rfr, std::vector<double>& gearing_1,
    std::vector<double>& gearing_2,
    const std::vector<const std::vector<double>>& volatilities,
    const double initial_val,
    const std::shared_ptr<mc::random_system_setting> rsconfig)
    : rfr{rfr},
      gearing_1{gearing_1},
      gearing_2{gearing_2},
      volatilities{volatilities},
      mc::random_process_update_base_setting(dim, initial_val, rsconfig){

      };

// Why default contructor of bs_update_setting is called..?
bs_random_process::bs_random_process(const bs_update_setting & setting){

};
}  // namespace qlib::bs::mc