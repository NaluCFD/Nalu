/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef PulseVelocityAuxFunction_h
#define PulseVelocityAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class PulseVelocityAuxFunction : public AuxFunction
{
public:

  PulseVelocityAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    const std::vector<double> &params);

  virtual ~PulseVelocityAuxFunction() {}
  
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
  double ux_;
  double uy_;
  double uz_;
  double pulseTime_;
};

} // namespace nalu
} // namespace Sierra

#endif
