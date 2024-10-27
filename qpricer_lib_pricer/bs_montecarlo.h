#pragma once
#include <memory>
#include <vector>

#include "../qpricer_lib_core/montecarlo.hpp"

namespace qlib::bs::mc {
	namespace mc = qlib::core::mc;

	struct bs_update_setting final : mc::random_process_update_base_setting {
		bs_update_setting() = delete;
		bs_update_setting(int dim, double rfr, std::vector<double>& gearing_1,
			std::vector<double>& gearing_2,
			const std::vector<double>& volatilities,
			const double initial_val,
			const std::shared_ptr<mc::random_system_setting> rsconfig);
		const std::vector<double> gearing_1;
		const std::vector<double> gearing_2;
		const double rfr;
		const std::vector<double> volatilities;
		const uint_fast64_t steps;
	};

	class bs_random_process : mc::i_random_process<bs_update_setting> {
	public:
		bs_random_process() = delete;
		bs_random_process(const bs_update_setting& setting);

		mc::path_t get_path(
			const std::vector<std::vector<double>>& uncertainties_increments,
			const std::vector<double>& time_increments, const int steps) override;

	private:
		void _update(const std::vector<double> uncertainties_increment,
			const bs_update_setting& setting, const double time_increment, const int current_step) override;
		bs_update_setting _setting;
		double current_val;
		mc::path_t path;
	};

	template <class ...BsRps>
	class bs_random_system : mc::i_random_system<BsRps ...> {
	public:
		bs_random_system() = delete;
		bs_random_system(const std::tuple<BsRps...>& rps, const mc::random_system_setting&);
	private:
		double payoff(std::tuple<BsRps...>& rps) override;
	};

}  // namespace qlib::bs::mc