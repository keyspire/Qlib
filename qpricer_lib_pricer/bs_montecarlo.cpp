#include "bs_montecarlo.h"
#include <numeric>
#include <cmath>


namespace qlib::bs::mc {

	bs_update_setting::bs_update_setting(
		int dim, double rfr, std::vector<double>& gearing_1,
		std::vector<double>& gearing_2,
		const std::vector<double>& volatilities,
		const double initial_val,
		const std::shared_ptr<mc::random_system_setting> rsconfig)
		: rfr{ rfr },
		gearing_1{ gearing_1 },
		gearing_2{ gearing_2 },
		volatilities{ volatilities },
		mc::random_process_update_base_setting(dim, initial_val, rsconfig), steps{ rsconfig->mcconfig->num_steps } {

	};



	bs_random_process::bs_random_process(const bs_update_setting& setting) : current_val{ setting.initial_val }, _setting{ setting } {
	};

	mc::path_t bs_random_process::get_path(const std::vector<std::vector<double>>& uncertainties_increments,
		const std::vector<double>& time_increments, const int steps) {
		path.resize(steps + 1);
		path[0] = current_val;
		for (uint_fast64_t i = 0; i < steps; ++i) {
			_update(uncertainties_increments[i], _setting, time_increments[i], i);
			path[i + 1] = current_val;
		}
		return path;
	}

	void bs_random_process::_update(const std::vector<double> uncertainties_increment,
		const bs_update_setting& setting, const double time_increment, const int current_step) {
		current_val = setting.rfr * setting.gearing_1.at(current_step) * time_increment + std::sqrt(time_increment) * std::inner_product(setting.volatilities.begin(), setting.volatilities.end(), uncertainties_increment.begin(), 0.);
		/**
		std::inner_product(setting.volatilities.begin(), setting.volatilities.end(), uncertainties_increment.begin(), 0.);*/

	}

	//double bs_random_system_1_w_mirror::payoff(std::tuple<bs_random_process, bs_random_process>& rps) {
	//		
	//}




}  // namespace qlib::bs::mc