/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/FixedHeightMixFracAuxFunction.h>
#include <algorithm>
#include <NaluEnv.h>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

FixedHeightMixFracAuxFunction::FixedHeightMixFracAuxFunction(
    const std::vector<double> &theParams) :
  AuxFunction(0,1),
  hY_(0.0),
  dY_(0.1),
  amp_(0.0),
  freq_(0.0)
{
  // extract the params - if they are supplied (optional)
  if ( theParams.size() > 0) 
    hY_ = theParams[0];
  if ( theParams.size() > 1)
    dY_ = theParams[1];
  if ( theParams.size() > 2)
    amp_ = theParams[2];
  if ( theParams.size() > 3)
    freq_ = theParams[3];

}

void
FixedHeightMixFracAuxFunction::do_evaluate(
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

    const double y = coords[1];
    const double x = coords[0];
    double dy; 
    if(std::abs(amp_)>0.0000001)
      dy = y-amp_*sin(freq_*x);
    else
      dy = y-hY_;
         
    fieldPtr[0] = 0.5*(std::erf(dy/dY_)+1.0);

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
