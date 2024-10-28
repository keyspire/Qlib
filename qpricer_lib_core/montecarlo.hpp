#pragma once
#include <memory>
#include <random>
#include <type_traits>
#include <vector>

#include "./template/template.h"

namespace qlib::core::mc {
using path_t = std::vector<double>;
using relation_vectors_t = std::vector<std::vector<uint_fast64_t>>;

struct mc_setting {
  mc_setting() = delete;
  mc_setting(const std::vector<double>& t_incs, const uint_fast64_t num_steps,
             uint_fast64_t num_runs);
  const std::vector<double> t_incs;
  const uint_fast64_t num_steps;
  const uint_fast64_t num_runs;
  std::random_device seed_gen;
};

struct random_system_setting {
  random_system_setting() = delete;
  random_system_setting(const relation_vectors_t& rel_vecs,
                        const uint_fast64_t num_uncertainties,
                        const std::shared_ptr<mc_setting>& mcconfig);
  // rel_vecs defines the number of unceratainty sources in random processes and
  // their relation. e.g., if there are 5 uncertainty sources (dW1, dW2, dW3,
  // dW4, dW5) and random processes X1, X2 depends on them as dX1 (dW1, dW3),
  // dX2 (dW2, dW4, dW5) we have rel_vecs = {{0, 2}, {1, 3, 4}}. NOTE: generally
  // rel_vecs is not a matrix.
  // TODO: implement factory method which ensures consistency all of the
  // settings.
  const relation_vectors_t rel_vecs;
  // The number of uncertainties included in the system.
  // In the above case, num_uncertainties = 5.
  const uint_fast64_t num_uncertainties;
  const std::shared_ptr<mc_setting> mcconfig;
};

struct random_process_update_base_setting {
  random_process_update_base_setting() = delete;
  random_process_update_base_setting(
      const double initial_val, const uint_fast64_t index_in_system,
      const std::shared_ptr<random_system_setting> rsconfig);
  // This pointer rsconfig ensures that a random process cannot be instantiated
  // without a wrapping random system.
  const std::shared_ptr<random_system_setting> rsconfig;
  // index in a random system.
  const uint_fast64_t index_in_system;
  const double initial_val;
};

class i_random_process {
 public:
  virtual ~i_random_process() = 0;

  virtual path_t get_path(
      const std::vector<std::vector<double>>& uncertainties_increments,
      const std::vector<double>& time_increments,
      const uint_fast64_t steps) = 0;

 private:
  virtual void _update(const std::vector<double>& uncertainties_increment,
                       const double time_increment,
                       const uint_fast64_t current_step) = 0;
};
// Ç±ÇÃÇ‹Ç‹Ç≈ÇÕ payoff Ç™é¿ëïÇ≈Ç´Ç»Ç¢ÅcÅB
template <class... Rps>
class i_random_system {
  static_assert(qlib::utils::core::all_base_of<i_random_process, Rps>(),
                "Random Process Setting(s) must be properly given.");

 public:
  virtual ~i_random_system() = 0;
  virtual std::vector<path_t> get_paths() = 0;
};

template <class Rs, class... Rps>
class i_mc {
  static_assert((std::is_base_of<i_random_system<Rps...>, Rs>::value),
                "Random Process Setting(s) must be properly given.");

 public:
  ~i_mc() = delete;
  /*mc(const i_random_system<>& ix);*/
  virtual void run(const mc_setting& setting) = 0;
  virtual double payoff(Rs& rs) = 0;
};

}  // namespace qlib::core::mc
