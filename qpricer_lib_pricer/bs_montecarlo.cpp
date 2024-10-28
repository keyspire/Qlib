#include "bs_montecarlo.h"

#include <cmath>
#include <numeric>
#include <random>

#include "random/random.h"

namespace qlib::bs::mc {

bs_update_setting::bs_update_setting(
    double rfr, const std::vector<double>& gearing_1, const std::vector<double>& gearing_2,
    const std::vector<double>& volatilities, const double initial_val,
    const std::shared_ptr<mc::random_system_setting> rsconfig,
    const uint_fast64_t index_in_system)
    : rfr{rfr},
      gearing_1{gearing_1},
      gearing_2{gearing_2},
      volatilities{volatilities},
      mc::random_process_update_base_setting(initial_val, index_in_system,
                                             rsconfig),
      steps{rsconfig->mcconfig->num_steps} {

      };

bs_random_process::bs_random_process(const bs_update_setting& setting)
    : current_val{setting.initial_val}, _setting{setting} {};

mc::path_t bs_random_process::get_path(
    const std::vector<std::vector<double>>& uncertainties_increments,
    const std::vector<double>& time_increments, const uint_fast64_t steps) {
  path.resize(steps + 1);
  path[0] = current_val;
  for (uint_fast64_t i = 0; i < steps; ++i) {
    _update(uncertainties_increments[i], time_increments[i], i);
    path[i + 1] = current_val;
  }
  return path;
}

void bs_random_process::_update(
    const std::vector<double>& uncertainty_increments,
    const double time_increment, const uint_fast64_t current_step) {
  double _uncertainty_part{0.};

  // a vector uncertainty_increments is a time-slice of
  // uncertainties_increments. (dW_1(t_1), dW_2(t_1), ..., dW_n(t_1)) =
  // time-slice at t_1 = uncertainty-increments. (time-slice at t_1, time-slice
  // at t_2, ..., ) = uncertain"ties" increments. Since the relation may specify
  // a discontinuous elements in the uncertainty_increments, inner-product
  // cannot be used below.
  for (uint_fast64_t i = 0, _size = _setting.volatilities.size(); i < _size;
       ++i) {
    _uncertainty_part +=
        _setting.volatilities.at(i) *
        uncertainty_increments.at(
            _setting.rsconfig->rel_vecs.at(_setting.index_in_system).at(i));
  }
  current_val +=
      _setting.rfr * _setting.gearing_1.at(current_step) * time_increment +
      std::sqrt(time_increment) * _setting.gearing_1.at(current_step) *
          _uncertainty_part;
  /* std::inner_product(_setting.volatilities.begin(),
                      _setting.volatilities.end(),
                      uncertainties_increment.begin(), 0.);*/
}

template <class... BsRps>
bs_random_system<BsRps...>::bs_random_system(
    const std::tuple<BsRps...>& rps, const mc::random_system_setting& setting)
    : _rps{rps}, _setting{setting} {
  //// erase unnecessary columns in the relation matrix
  // for (uint_fast64_t j = 0, _column_size = _rel_mat_trimmed.at(0).size(),
  //                    _row_size = _rel_mat_trimmed.size();
  //      j < _column_size; ++j) {
  //   uint_fast64_t _sum = 0;
  //   for (uint_fast64_t i = 0; i < _row_size; ++i) {
  //     _sum += _rel_mat_trimmed.at(i).at(j);
  //   }
  //   if (_sum == 0) {
  //     for (uint_fast64_t i = 0; i < _row_size; ++i) {
  //       _rel_mat_trimmed.at(i).erase(_rel_mat_trimmed.at(i).begin() + j);
  //     }
  //   }
  // }
}

// template <class... BsRps>
// template <int N>
// mc::path_t bs_random_system<BsRps...>::get_ith_path() {
//   return std::get<N>(_rps).get_path();
// }

template <class... BsRps>
std::vector<mc::path_t> bs_random_system<BsRps...>::get_paths() {
  // init seed
  _seed = _setting.mcconfig->seed_gen();
  std::mt19937_64 _engine = qlib::utils::core::create_engine(_seed);

  // TODO
  for (uint_fast64_t i = 0, _num_steps = _setting.mcconfig->num_steps;
       i < _setting.num_uncertainties; ++i) {
    _uncertainties_increments.at(i).resize(_num_steps);
    for (uint_fast64_t j = 0; j < _num_steps; ++j) {
      _uncertainties_increments.at(i).push_back(_dist(_engine));
    }
  }

  // args are assumed as rps
  std::apply(
      [this](auto&&... args) {
        (_paths.push_back(args.get_path(_uncertainties_increments,
                                        _setting.mcconfig->t_incs,
                                        _setting.mcconfig->num_steps)),
         ...);
      },
      _rps);

  // another implementation

  return _paths;
}

//template <class... BsRps>
//template <int N>
//void bs_random_system<BsRps...>::_register_path_for_each_tuples() {
//  _path.push_back()
//
//      return;
//}

}  // namespace qlib::bs::mc