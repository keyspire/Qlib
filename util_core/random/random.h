#pragma once
#include <random>

namespace qlib::utils::core {
	inline std::mt19937_64 create_engine(uint_fast64_t seed) {
		std::mt19937_64 _engine(seed);
		return _engine;
	}
}