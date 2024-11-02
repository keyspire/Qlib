#pragma once
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

namespace analyze {
using result_t = std::vector<double>;

template <typename T>
std::vector<std::vector<T>> get_partition(const std::vector<T>& vec,
                                          const uint64_t cell_size) {
  const size_t length = vec.size();
  auto size = length / cell_size;
  auto is_zero_remainder = (length % cell_size == 0) ? 0 : 1;
 std::vector<std::vector<T>> div_vec;
  //div_vec.resize(size + is_zero_remainder);
  for (auto i = 0; i < size + is_zero_remainder; ++i) {
    auto _temp_vec = std::vector<T>{};
    if (i == size) {
      for (auto j = 0; j < length % cell_size; ++j) {
        _temp_vec.push_back(vec.at(cell_size * i + j));
      }
    } else if (i < size) {
      for (auto j = 0; j < cell_size; ++j) {
        _temp_vec.push_back(vec.at(cell_size * i + j));
      }
    }
    div_vec.push_back(_temp_vec);
  }
  return div_vec;
}


template <typename T>
T variance(const std::vector<T>& vec) {
  const size_t sz = vec.size();
  if (sz <= 1) {
    return 0.0;
  }

  // Calculate the mean
  const T mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

  // Now calculate the variance
  auto variance_func = [&mean, &sz](T accumulator, const T& val) {
    return accumulator + ((val - mean) * (val - mean) / (sz - 1));
  };

  return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

std::vector<result_t> get_95p_confidence_interval_t_test(const std::vector<result_t> divided_results);

class feed_normal_rv {
 public:
  feed_normal_rv(const uint64_t seed);
  double get_single_normal_rv();

 private:
  std::mt19937_64 _engine;
  std::normal_distribution<> _dist;
};

}  // namespace main
