#pragma once
#include "bs_montecarlo.h"
namespace qlib::bs::mc {
mc::mc_setting mcconfig = mc::mc_setting{{1.0, 1.0}, 2, 100};

mc::random_system_setting rs_setting = mc::random_system_setting{
    {{1}, {1}}, 1, std::make_shared<mc::mc_setting>(mcconfig)};

bs_update_setting setting1 =
    bs_update_setting{0.1,
                      std::vector<double>{1.0, 1.0},
                      std::vector<double>{1.0, 1.0},
                      std::vector<double>{0.1},
                      0.0,
                      std::make_shared<mc::random_system_setting>(rs_setting),
                      0};

bs_update_setting setting2 =
    bs_update_setting{0.1,
                      std::vector<double>{1.0, 1.0},
                      std::vector<double>{-1.0, -1.0},
                      std::vector<double>{0.1},
                      0.0,
                      std::make_shared<mc::random_system_setting>(rs_setting),
                      0};

bs_random_system<bs_random_process, bs_random_process> rs_fx_option{
    std::tuple<bs_random_process, bs_random_process>{
        bs_random_process{setting1}, bs_random_process{setting2}},
    mc::random_system_setting{rs_setting}};

template <class... BsRps>
class mc_fx_option : mc::i_mc {
  static_assert(
      qlib::utils::core::all_base_of<bs_random_process, BsRps...>(),
      "mc_fx_option is currently a class whose random processes are only "
      "of Black-Scholes type.");

 public:
  mc_fx_option() = delete;
  mc_fx_option(mc::mc_setting mcconfig, bs_random_system<BsRps...>, std::string flag);
  void run() override;

 private:
  double payoff() override;

 private:
  mc::mc_setting _setting;
  bs_random_system<BsRps...> _rs;
};
}  // namespace qlib::bs::mc