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
};

// T must be derived from random_process_update_base_setting
template <class T, typename = std::enable_if_t<std::is_base_of<
                       random_process_update_base_setting, T>::value>>
class i_random_process {
 public:
  virtual ~i_random_process() = 0;

  virtual path_t get_path(
      const std::vector<const std::vector<double>>& uncertainties_increments,
      const std::vector<double>& time_increments) = 0;

 private:
  virtual void _update(const std::vector<double> uncertainties_increment,
                       const T& setting) = 0;
};
// このままでは payoff が実装できない…。
template <class... SettingRpUpdate>
class i_random_system {
  static_assert((std::is_base_of<random_process_update_base_setting,
                                 SettingRpUpdate>::value &&
                 ...),
                "Random Process Setting(s) must be properly given.");

 public:
  virtual ~i_random_system() = 0;
  // vector の中身ポインタにしないと object slicing ?
  virtual double payoff(
      std::vector<i_random_process<SettingRpUpdate>...>& rps) = 0;

 private:
  virtual std::vector<path_t> get_paths() = 0;

  virtual path_t get_ith_path(const int i) = 0;
};

template <class SettingRpUpdate>
class mc {
 public:
  mc() = delete;
  std::vector<double> run(
      const mc_setting&,
      i_random_system<random_system_setting, SettingRpUpdate>& rs);
};
}  // namespace qlib::core::mc
