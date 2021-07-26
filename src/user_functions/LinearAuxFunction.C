/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/LinearAuxFunction.h>
#include <MeshMotionInfo.h>
#include <PecletFunction.h>
#include <Realm.h>
#include <SolutionOptions.h>

// basic c++
#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

LinearAuxFunction::LinearAuxFunction(
  const std::vector<double> &params,
  const unsigned beginPos,
  const unsigned endPos) :
  AuxFunction(beginPos, endPos),
  m_(0.0),
  b_(0.0),
  coordPd_(0),
  fieldPd_(0)
{
  // throw is not fully specified
  if ( params.size() != 6 )
    throw std::runtime_error("LinearAuxFunction::error: Requires six parameters");
  
  const double qL = params[0];
  const double qR = params[1];
  const double xL = params[2];
  const double xR = params[3];
  coordPd_ = params[4];
  fieldPd_ = params[5];

  // compute slope and intercept
  m_ = (qR - qL)/(xR - xL);
  b_ = qL - m_*xL;
}

void
LinearAuxFunction::do_evaluate(
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
    
    // set all to zero
    for ( unsigned j = 0; j < fieldSize; ++j )
      fieldPtr[j] = 0.0;

    // populate for principle direction
    const double y = m_*coords[coordPd_] + b_;
    fieldPtr[fieldPd_] = y;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
