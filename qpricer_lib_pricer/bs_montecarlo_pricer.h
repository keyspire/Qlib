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

}  // namespace qlib::bs::mc