#pragma once
#include <type_traits>

namespace qlib::utils::core{
	template<typename Base, typename... Args>
	constexpr auto all_base_of()
	{
		return (std::is_base_of<Base, Args>::value && ...);
	}
}