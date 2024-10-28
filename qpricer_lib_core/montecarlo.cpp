#include "montecarlo.hpp"

#include <exception>
namespace qlib::core::mc {
random_process_update_base_setting::random_process_update_base_setting(
    const double initial_val, const uint_fast64_t index_in_system,
    const std::shared_ptr<random_system_setting> rsconfig)
    : rsconfig{rsconfig},
      initial_val{initial_val},
      index_in_system{index_in_system} {};

mc_setting::mc_setting(const std::vector<double>& t_incs,
                       const uint_fast64_t num_steps,
                       const uint_fast64_t num_runs)
    : t_incs{t_incs},
      num_steps{num_steps},
      num_runs{num_runs} {

      };

random_system_setting::random_system_setting(
    const relation_vectors_t& rel_vecs, const uint_fast64_t num_uncertainties,
    const std::shared_ptr<mc_setting>& mcconfig)
    : rel_vecs{rel_vecs},
      num_uncertainties{num_uncertainties},
      mcconfig{mcconfig} {
  for (uint_fast64_t i = 0, _row_size = rel_vecs.size(); i < _row_size; ++i) {
    for (uint_fast64_t j = 0, _column_size = rel_vecs.at(i).size();
         i < _column_size; ++j) {
      if (rel_vecs.at(i).at(j) >= num_uncertainties) {
        throw std::logic_error(
            "The relation vectors and the number of uncertainties are "
            "inconsistent.");
      }
    }
  }
}
}  // namespace qlib::core::mc