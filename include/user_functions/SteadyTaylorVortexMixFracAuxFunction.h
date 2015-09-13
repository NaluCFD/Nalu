/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SteadyTaylorVortexMixFracAuxFunction_h
#define SteadyTaylorVortexMixFracAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class SteadyTaylorVortexMixFracAuxFunction : public AuxFunction
{
public:

  SteadyTaylorVortexMixFracAuxFunction();

  virtual ~SteadyTaylorVortexMixFracAuxFunction() {}
  
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
  const double znot_;
  const double amf_;
  const double pi_;
};

} // namespace nalu
} // namespace Sierra

#endif
