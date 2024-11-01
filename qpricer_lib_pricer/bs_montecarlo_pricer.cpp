#include "bs_montecarlo_pricer.h"

#include <exception>
#include "bs_montecarlo.h"
namespace qlib::bs::mc {
template <class... BsRps>
double mc_fx_option<BsRps...>::payoff() {
  throw std::logic_error("payoff() is not yet implemented");
}

template <>
double mc_fx_option<bs_random_process, bs_random_process>::payoff() {
  if (_setting.flag == "negative_corr") {
    const auto _res = _rs.get_paths();
    return (_res->at(0).back(), 0);//TODO
  } else {
    throw std::logic_error("payoff() is not yet implemented");
  }
}
}  // namespace qlib::bs::mc