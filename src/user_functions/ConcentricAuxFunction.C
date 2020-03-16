/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/ConcentricAuxFunction.h>
#include <algorithm>
#include <NaluEnv.h>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

ConcentricAuxFunction::ConcentricAuxFunction() :
  AuxFunction(0,1),
  ti_(300.0),
  //t1_(302.834),
  t1_(302.83443018241820254843332804739),
  ri_(0.03),
  r1_(0.05)
{
  // nothing
}

void
ConcentricAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];
    
    const double R = std::sqrt(x*x + y*y + z*z);
    const double temp = ti_ + (t1_-ti_)*(1/ri_ - 1/R)/(1.0/ri_ - 1/r1_);
    fieldPtr[0] = temp;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
