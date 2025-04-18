/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef RecMixFracAuxFunction_h
#define RecMixFracAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class RecMixFracAuxFunction : public AuxFunction
{
public:

  RecMixFracAuxFunction();

  virtual ~RecMixFracAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;  
};

} // namespace nalu
} // namespace Sierra

#endif
