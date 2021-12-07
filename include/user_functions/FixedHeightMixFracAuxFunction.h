/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef FixedHeightMixFracAuxFunction_h
#define FixedHeightMixFracAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class FixedHeightMixFracAuxFunction : public AuxFunction
{
public:

  FixedHeightMixFracAuxFunction(
    const std::vector<double> &theParams);

  virtual ~FixedHeightMixFracAuxFunction() {}
  
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
  double hY_;
  double dY_;
};

} // namespace nalu
} // namespace Sierra

#endif
