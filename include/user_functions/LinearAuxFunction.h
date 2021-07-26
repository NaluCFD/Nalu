/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef LinearAuxFunction_h
#define LinearAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class LinearAuxFunction : public AuxFunction
{
public:

  LinearAuxFunction(
    const std::vector<double> &params,
    const unsigned beginPos,
    const unsigned endPos);

  virtual ~LinearAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
private:
  double m_;
  double b_;
  int coordPd_;
  int fieldPd_;
};

} // namespace nalu
} // namespace Sierra

#endif
