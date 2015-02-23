#include <vector>
#include <cstddef>
#include <algorithm>

#include "fitter.cc"

namespace Numina {

class Trace {
  public:
  Trace() {}
  void push_back(double x, double y) {
    xtrace.push_back(x);
    ytrace.push_back(y);

  }

  std::vector<double> xtrace;
  std::vector<double> ytrace;

  void reverse() {
    std::reverse(xtrace.begin(), xtrace.end());
    std::reverse(ytrace.begin(), ytrace.end());
  }

  double predict(double x) const {

    size_t n = std::min<size_t>(5, xtrace.size());
    LinearFit mm = linear_fitter(xtrace.end() - n, xtrace.end(), ytrace.end() - n, ytrace.end());
    return mm.slope * x + mm.intercept;
  }
};

} // namespace numina
