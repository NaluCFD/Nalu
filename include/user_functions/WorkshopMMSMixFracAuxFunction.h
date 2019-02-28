/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef WorkshopMMSMixFracAuxFunction_h
#define WorkshopMMSMixFracAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class WorkshopMMSMixFracAuxFunction : public AuxFunction
{
public:

  WorkshopMMSMixFracAuxFunction();

  virtual ~WorkshopMMSMixFracAuxFunction() {}
  
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
  const double u_;
  const double v_;
  const double sigma_;
};

} // namespace nalu
} // namespace Sierra

#endif
