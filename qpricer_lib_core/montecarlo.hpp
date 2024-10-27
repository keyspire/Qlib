#pragma once
#include <memory>
#include <type_traits>
#include <vector>

#include "../qpricer_lib/pch.h"

namespace qlib::core::mc {
	using path_t = std::vector<double>;
	using relation_matrix_t = std::vector<std::vector<double>>;

	struct mc_setting {
		mc_setting() = delete;
		const std::vector<double> t_incs;
		const uint_fast64_t num_steps;
		const uint_fast64_t num_runs;
	};

	struct random_system_setting {
		random_system_setting() = delete;
		random_system_setting(const relation_matrix_t& rel_mat,
			const std::shared_ptr<mc_setting> mcconfig);
		const relation_matrix_t rel_mat;
		const std::shared_ptr<mc_setting> mcconfig;
	};

	struct random_process_update_base_setting {
		random_process_update_base_setting() = delete;
		random_process_update_base_setting(
			const int dim, const double initial_val,
			const std::shared_ptr<random_system_setting> rsconfig);
		// number of uncertainty sources
		const int dim;
		const std::shared_ptr<random_system_setting> rsconfig;
		const double initial_val;
		const std::vector<double> t_incs;
	};

	// T must be derived from random_process_update_base_setting
	template <class T, typename = std::enable_if_t<std::is_base_of<
		random_process_update_base_setting, T>::value>>
		class i_random_process {
		public:
			virtual ~i_random_process() = 0;

			virtual path_t get_path(
				const std::vector<std::vector<double>>& uncertainties_increments,
				const std::vector<double>& time_increments,
				const int steps) = 0;

		private:
			virtual void _update(std::vector<double> uncertainties_increment,
				const T& setting, const double time_increment, const int current_step) = 0;
	};
	// このままでは payoff が実装できない…。
	template <class... Rps>
	class i_random_system {
		static_assert((std::is_base_of<i_random_process,
			Rps>::value &&
			...),
			"Random Process Setting(s) must be properly given.");

	public:
		virtual ~i_random_system() = 0;
		// vector の中身ポインタにしないと object slicing ?

	private:
		virtual std::vector<path_t> get_paths() = 0;

		virtual path_t get_ith_path(const int i) = 0;
	};

	template <class Rs, class... Rps>
	class mc {
		static_assert((std::is_base_of<i_random_system<Rps...>, Rs>::value),
			"Random Process Setting(s) must be properly given.");
	public:
		mc() = delete;
		/*mc(const i_random_system<>& ix);*/
		void run();
		virtual double payoff(
			Rs& rs) = 0;
	};

}  // namespace qlib::core::mc
