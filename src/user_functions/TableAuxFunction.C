/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/TableAuxFunction.h>
#include <MeshMotionInfo.h>
#include <PecletFunction.h>
#include <SolutionOptions.h>

// basic c++
#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

TableAuxFunction::TableAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  std::vector<double> params) :
  AuxFunction(beginPos, endPos),
  fieldComp_(0),
  interpComp_(0),
  tableOffset_(2)
{
  /* sample: 0, 1, 
     0.0, 0.0, --> first two entries taken
     1.0, 2.0  --> begin table data; in pairs 
  */
  fieldComp_ = params[0];
  interpComp_ = params[1];
  for ( size_t k = tableOffset_; k < params.size(); k += 2 ) {
    table_.push_back(std::make_pair(params[k+0],params[k+1]));
  }

  // sort the table - just in case
  std::sort(table_.begin(), table_.end());
  
  // make sure it's more than one entry
  if ( table_.size() < 2 )
    throw std::runtime_error("TableAuxFunction::Error() Table must hold at least two entries");
}

TableAuxFunction::~TableAuxFunction()
{
  // nothing
}

void
TableAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    // create the pair, coordinates or time
    const std::pair<double,double> cPair = (interpComp_ < 3 )
      ? std::make_pair(coords[interpComp_],0.0) 
      : std::make_pair(time,0.0);

    // interpolate
    const double interpValue = interp(cPair);
    
    // first assign to zero
    for ( unsigned i = 0; i < fieldSize; ++i )
      fieldPtr[i] = 0.0;
    
    // assign
    fieldPtr[fieldComp_] = interpValue;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

double
TableAuxFunction::interp(const std::pair<double,double> &x) const
{
  double returnValue = 0.0;
  if ( x.first <= table_[0].first ) {
    returnValue = table_[0].second;
  }
  else if ( x.first >= table_[table_.size()-1].first ) {
    returnValue = table_[table_.size()-1].second;
  }
  else {
    // find the lower bound
    std::vector<std::pair<double, double> >::const_iterator it = std::lower_bound( table_.begin(), table_.end(), x );    
    if ( it == table_.end() ) {
      throw std::runtime_error("TableAuxFunction::Error() Table entry not found....");
    }
    else {
      const double dx = (*it).first - (*(it-1)).first;
      const double dy = (*it).second - (*(it-1)).second;
      returnValue = (*(it-1)).second + (x.first  - (*(it-1)).first) * dy / dx;
    }
  }
  return returnValue;
}

} // namespace nalu
} // namespace Sierra
