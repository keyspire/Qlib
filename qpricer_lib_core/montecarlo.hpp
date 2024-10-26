#pragma once
#include <type_traits>
#include <vector>

#include "../qpricer_lib/pch.h"

namespace qlib::core::mc {
using path_t = std::vector<double>;
using relation_matrix_t = std::vector<std::vector<double>&>;

struct mc_setting {
  mc_setting() = delete;
  const std::vector<double> t_incs;
  const uint_fast64_t num_steps;
  const uint_fast64_t num_runs;
};

struct random_system_base_setting : mc_setting {
  random_system_base_setting() = delete;
  const relation_matrix_t rel_mat;
};

struct random_process_update_base_setting : random_system_base_setting {
  random_process_update_base_setting() = delete;
  // number of uncertainty sources
  const int dim;
};

// T must be derived from random_process_update_base_setting
template <class T, typename = std::enable_if_t<std::is_base_of<
                          random_process_update_base_setting, T>::value>>
class i_random_process {
 public:
  ~i_random_process() = delete;

  virtual path_t get_path(
      const std::vector<const std::vector<double>&>& uncertainties_increments,
      const std::vector<double>& time_increments) = 0;

 private:
  virtual void _update(const std::vector<double> uncertainties_increment,
                       const T& setting) = 0;
};

template <class SettingRs, class SettingRpUpdate,
          typename = std::enable_if<
              std::is_base_of<random_system_base_setting, SettingRs>::value &&
              std::is_base_of<random_process_update_base_setting,
                              SettingRpUpdate>::value>>
class i_random_system {
 public:
  ~i_random_system() = delete;
  //std::vector < i_random_process<SettingRpUpdate>&>
  virtual double payoff(std::vector<i_random_process<SettingRpUpdate>&>& rps,
                        const SettingRs& setting) = 0;

 private:
  virtual std::vector<path_t> get_paths() = 0;

  virtual path_t get_ith_path(const int i) = 0;
};

class i_mc {
 public:
  ~i_mc() = delete;
  virtual std::vector<double> run(mc_setting&) = 0;
};
}  // namespace qlib::core::mc
