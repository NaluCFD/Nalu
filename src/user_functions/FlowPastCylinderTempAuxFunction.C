/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/FlowPastCylinderTempAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

FlowPastCylinderTempAuxFunction::FlowPastCylinderTempAuxFunction() :
  AuxFunction(0,1),
  h_(0.420),
  k_(1.42690e-2),
  pi_(acos(-1.0)),
  iMin_(0),
  iMax_(24)
{
  // initialize the data; Angle (degree) and Temperature (C) 
  experimentalData_[0][0]  =  0.00000000000; experimentalData_[0][1]  = 34.8567771912;
  experimentalData_[1][0]  =  15.0000000000; experimentalData_[1][1]  = 34.9078674316;
  experimentalData_[2][0]  =  30.0000000000; experimentalData_[2][1]  = 35.6716117859;
  experimentalData_[3][0]  =  45.0000000000; experimentalData_[3][1]  = 36.5256385803;
  experimentalData_[4][0]  =  60.0000000000; experimentalData_[4][1]  = 37.9703292847;
  experimentalData_[5][0]  =  75.0000000000; experimentalData_[5][1]  = 38.8416900635;
  experimentalData_[6][0]  =  90.0000000000; experimentalData_[6][1]  = 40.3171539307;
  experimentalData_[7][0]  = 105.0000000000; experimentalData_[7][1]  = 42.0815200806;
  experimentalData_[8][0]  = 120.0000000000; experimentalData_[8][1]  = 41.2935409546;
  experimentalData_[9][0]  = 135.0000000000; experimentalData_[9][1]  = 40.9183769226;
  experimentalData_[10][0] = 150.0000000000; experimentalData_[10][1] = 40.2051162720;
  experimentalData_[11][0] = 165.0000000000; experimentalData_[11][1] = 40.3736419678;
  experimentalData_[12][0] = 180.0000000000; experimentalData_[12][1] = 40.3435249329;
  experimentalData_[13][0] = 195.0000000000; experimentalData_[13][1] = 40.5962333679;
  experimentalData_[14][0] = 210.0000000000; experimentalData_[14][1] = 40.7046546936;
  experimentalData_[15][0] = 225.0000000000; experimentalData_[15][1] = 40.5719375610;
  experimentalData_[16][0] = 240.0000000000; experimentalData_[16][1] = 39.5904731750;
  experimentalData_[17][0] = 255.0000000000; experimentalData_[17][1] = 38.0906562805;
  experimentalData_[18][0] = 270.0000000000; experimentalData_[18][1] = 36.7716522217;
  experimentalData_[19][0] = 285.0000000000; experimentalData_[19][1] = 35.9330558777;
  experimentalData_[20][0] = 300.0000000000; experimentalData_[20][1] = 35.2598419189;
  experimentalData_[21][0] = 315.0000000000; experimentalData_[21][1] = 34.9156837463;
  experimentalData_[22][0] = 330.0000000000; experimentalData_[22][1] = 34.4239540100;
  experimentalData_[23][0] = 345.0000000000; experimentalData_[23][1] = 34.2773284912;
  experimentalData_[24][0] = 360.0000000000; experimentalData_[24][1] = 34.8567771912;
}

void
FlowPastCylinderTempAuxFunction::do_evaluate(
  const double *coords,
  const double t,
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

    // determine the quadrant
    double addMe = 0.0;
    if ( x < h_ )
      if ( y > k_ )
        addMe = 0.0;   // Q1
      else
        addMe = 270.0; // Q3
    else
      if ( y > k_ )
        addMe = 90.0;  // Q2
      else
        addMe = 180.0; // Q4

    const double radius = std::sqrt((x-h_)*(x-h_) + (y-k_)*(y-k_));
    const double theta = 90.0/pi_*acos((x-h_)/radius) + addMe;
    
    // find the value via linear interpolation
    const double findT = interpolate_data(theta);
    
    // convert to K
    fieldPtr[0] = findT+273.15;
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

int
FlowPastCylinderTempAuxFunction::find_index( 
  const double z,
  int iMin,
  int iMax) const
{  
  // mid point of table
  int iMid = (iMin + iMax)/2;

  // find the mid point value for z; call it the candidate
  const double candidateZ = experimentalData_[iMid][0];
  
  if ( iMid < iMin) {
    return iMid;
  }

  if ( candidateZ > z ) {    
    // in the lower
    return find_index(z, iMin, iMid-1);
  }
  else {
    // in the upper
    return find_index(z, iMid+1, iMax);
  }
} 

double
FlowPastCylinderTempAuxFunction::interpolate_data( 
  const double z) const
{
  // find the index
  int foundIndex = find_index(z, iMin_, iMax_);

  // proceed with the simple linear interpolation
  double interpProp = 0.0;

  // check for high/low bounds
  if ( foundIndex == iMin_ || foundIndex == iMax_)
    interpProp = experimentalData_[foundIndex][1];
  else {
    // not on the outer edge of the table
    if ( experimentalData_[foundIndex][0] > z ) {
      // on the high side
      int foundM1 = foundIndex-1;
      
      // values of mixFrac
      const double x1 = experimentalData_[foundIndex][0];
      const double x0 = experimentalData_[foundM1][0];
      
      // property
      const double y1 = experimentalData_[foundIndex][1];
      const double y0 = experimentalData_[foundM1][1];
      interpProp = y0 + (y1 - y0 )*(z-x0)/(x1-x0);
    }
    else {
      // on the low side
      const int foundP1 = foundIndex+1;
      
      // values of mixFrac
      const double x0 = experimentalData_[foundIndex][0];
      const double x1 = experimentalData_[foundP1][0];
      
      // property
      const double y0 = experimentalData_[foundIndex][1];
      const double y1 = experimentalData_[foundP1][1];
      
      interpProp = y0 + (y1 - y0 )*(z-x0)/(x1-x0);
    }
  }
  return interpProp;
}

} // namespace nalu
} // namespace Sierra
