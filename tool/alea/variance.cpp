#include <mcanalyze_tools.hpp>

template <class X> template <class TimeseriesType>
typename return_type<TimeseriesType>::type impl_calculation<X>::calculate (TimeseriesType& timeseries) {
  return alps::alea::variance(timeseries);
}

impl_calculation<> calculation("Variance", "variance/value");

#include <mcanalyze_tools.ipp> 











