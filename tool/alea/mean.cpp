#include <mcanalyze_tools.hpp>

template <class X> template <class TimeseriesType>
typename return_type<TimeseriesType>::type impl_calculation<X>::calculate (TimeseriesType& timeseries) {
  return alps::alea::mean(timeseries);
}

impl_calculation<> calculation("Mean", "mean/value");

#include <mcanalyze_tools.ipp> 











