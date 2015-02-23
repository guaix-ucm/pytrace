#ifndef NU_FITTER_H
#define NU_FITTER_H
namespace Numina {

struct LinearFit {
  double slope;
  double intercept;
  bool fit;
};

template <class Iterator>
LinearFit linear_fitter(Iterator x1, Iterator x2, Iterator y1, Iterator y2);

} // namespace Numina

#endif // NU_FITTER_H
