#pragma once
#include <random>

class feed_normal_rv {
 public:
  feed_normal_rv(const uint64_t seed);
  double get_single_normal_rv();
 private:
  std::mt19937_64 _engine;
  std::normal_distribution<> _dist;
};
