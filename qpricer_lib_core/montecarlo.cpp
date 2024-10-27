#include "montecarlo.hpp"

namespace qlib::core::mc {
	random_process_update_base_setting::random_process_update_base_setting(
		const int dim, const double initial_val,
		const std::shared_ptr<random_system_setting> rsconfig)
		: dim{ dim }, rsconfig{ rsconfig }, initial_val{ initial_val } {};
}