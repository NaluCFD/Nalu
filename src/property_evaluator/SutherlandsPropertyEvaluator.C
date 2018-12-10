/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "property_evaluator/PropertyEvaluator.h"
#include "property_evaluator/SutherlandsPropertyEvaluator.h"
#include "property_evaluator/ReferencePropertyData.h"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <cmath>
#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// SutherlandsPropertyEvaluator - base class
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsPropertyEvaluator::SutherlandsPropertyEvaluator(
  const std::map<std::string, std::vector<double> > &polynomialCoeffsMap)
  : PropertyEvaluator()
{
  // extract the reference property data size
  size_t polySize = polynomialCoeffsMap.size();

  // resize polynomial coeffs
  polynomialCoeffs_.resize(polySize);

  // save off polynomial coeffs
  size_t k = 0;
  std::map<std::string, std::vector<double> >::const_iterator itpc;
  for ( itpc = polynomialCoeffsMap.begin();
        itpc!= polynomialCoeffsMap.end(); ++itpc, ++k) {
    std::vector<double> polyVec= (*itpc).second;
    size_t polyVecSize = polyVec.size();
    if ( polyVecSize < 3)
      throw std::runtime_error("Sutherlands polynomial evaluator needs three coeffs:");
    polynomialCoeffs_[k].resize(polyVecSize);
    double *pt_poly = &polynomialCoeffs_[k][0];
    for ( size_t j = 0; j < polyVecSize; ++j ) {
      pt_poly[j] = polyVec[j];
    }
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsPropertyEvaluator::~SutherlandsPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- set_reference_property_data -------------------------------------
//--------------------------------------------------------------------------
void
SutherlandsPropertyEvaluator::set_reference_property_data(
  const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap)
{
  // extract the reference property data size
  size_t refPropSize = referencePropertyDataMap.size();
  refMassFraction_.resize(refPropSize);

  // sanity check to ensure that reference and polynomialCoeffs size matches
  if ( refMassFraction_.size() != polynomialCoeffs_.size() )
    throw std::runtime_error("Sutherlands refMassFraction_ and polynomialCoeffs_ size do not match:");
    
  // save off reference values for Yk
  size_t k = 0;
  std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
  for ( itrp = referencePropertyDataMap.begin();
        itrp!= referencePropertyDataMap.end(); ++itrp, ++k) {
    ReferencePropertyData *propData = (*itrp).second;
    refMassFraction_[k] = propData->massFraction_;
  }
}

//--------------------------------------------------------------------------
//-------- compute_viscosity -----------------------------------------------
//--------------------------------------------------------------------------
double
SutherlandsPropertyEvaluator::compute_viscosity(
  const double &T,
  const double *pt_poly)
{
  const double muRef = pt_poly[0];
  const double TRef = pt_poly[1];
  const double SRef = pt_poly[2];
  
  return muRef*std::pow(T/TRef, 1.5)*(TRef+SRef)/(T+SRef);
}

//==========================================================================
// Class Definition
//==========================================================================
// SutherlandsYkrefPropertyEvaluator - evaluates mu based on temperature; Ykref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsYkrefPropertyEvaluator::SutherlandsYkrefPropertyEvaluator(
  const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
  const std::map<std::string, std::vector<double> > &polynomialCoeffsMap)
  : SutherlandsPropertyEvaluator(polynomialCoeffsMap)
{
  // set reference data
  set_reference_property_data(referencePropertyDataMap);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsYkrefPropertyEvaluator::~SutherlandsYkrefPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SutherlandsYkrefPropertyEvaluator::execute(
  double *indVarList,
  stk::mesh::Entity /*node*/)
{
  const double T = indVarList[0];
  
  // extract size and process sum
  const size_t ykSize = refMassFraction_.size();
  double sum_mu = 0.0;
  for ( size_t k = 0; k < ykSize; ++k ) {
    sum_mu += refMassFraction_[k]*compute_viscosity(T, &polynomialCoeffs_[k][0]);
  }

  return sum_mu;
}

//==========================================================================
// Class Definition
//==========================================================================
// SutherlandsYkPropertyEvaluator - evaluates mu based on temperature; Yk
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsYkPropertyEvaluator::SutherlandsYkPropertyEvaluator(
  const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
  stk::mesh::MetaData &metaData)
  : SutherlandsPropertyEvaluator(polynomialCoeffsMap),
    massFraction_(NULL),
    ykVecSize_(0)
{
  // sizing
  ykVecSize_ = polynomialCoeffsMap.size();
  
  // sanity check on sizes
  if ( polynomialCoeffs_.size() != ykVecSize_ )
    throw std::runtime_error("Sutherlands polynomialCoeffs_ and ykVecSize_ do not match:");
  
  // save off mass fraction field
  massFraction_ = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "mass_fraction");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsYkPropertyEvaluator::~SutherlandsYkPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SutherlandsYkPropertyEvaluator::execute(
  double *indVarList,
  stk::mesh::Entity node)
{
  const double T = indVarList[0];
  const double *massFraction = stk::mesh::field_data(*massFraction_, node);

  // process sum
  double sum_mu = 0.0;
  for ( size_t k = 0; k < ykVecSize_; ++k ) {
    sum_mu += massFraction[k]*compute_viscosity(T, &polynomialCoeffs_[k][0]);
  }

  return sum_mu;
}

//==========================================================================
// Class Definition
//==========================================================================
// SutherlandsYkTrefPropertyEvaluator - evaluates mu based on Yk and Tref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsYkTrefPropertyEvaluator::SutherlandsYkTrefPropertyEvaluator(
  const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
  stk::mesh::MetaData &metaData,
  const double tRef)
  : SutherlandsPropertyEvaluator(polynomialCoeffsMap),
    tRef_(tRef)
{
  // sizing
  ykVecSize_ = polynomialCoeffsMap.size();

  // sanity check on sizes
  if ( polynomialCoeffs_.size() != ykVecSize_ )
    throw std::runtime_error("Sutherlands polynomialCoeffs_ and ykVecSize_ do not match:");

  // save off mass fraction field
  massFraction_ = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "mass_fraction");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsYkTrefPropertyEvaluator::~SutherlandsYkTrefPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SutherlandsYkTrefPropertyEvaluator::execute(
  double */*indVarList*/,
  stk::mesh::Entity node)
{
  const double *massFraction = stk::mesh::field_data(*massFraction_, node);
  
  // process sum
  double sum_mu = 0.0;
  for ( size_t k = 0; k < ykVecSize_; ++k ) {
    sum_mu += massFraction[k]*compute_viscosity(tRef_, &polynomialCoeffs_[k][0]);
  }

  return sum_mu;
}

//==========================================================================
// Class Definition
//==========================================================================
// SutherlandsYkrefTrefPropertyEvaluator - evaluates mu based on reference
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsYkrefTrefPropertyEvaluator::SutherlandsYkrefTrefPropertyEvaluator(
  const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
  const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
  const double tRef)
  : SutherlandsPropertyEvaluator(polynomialCoeffsMap),
    tRef_(tRef)
{
  // set reference data
  set_reference_property_data(referencePropertyDataMap);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SutherlandsYkrefTrefPropertyEvaluator::~SutherlandsYkrefTrefPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
SutherlandsYkrefTrefPropertyEvaluator::execute(
  double */*indVarList*/,
  stk::mesh::Entity /*node*/)
{
  // extract size and process sum
  const size_t ykSize = refMassFraction_.size();
  double sum_mu = 0.0;
  for ( size_t k = 0; k < ykSize; ++k ) {
    sum_mu += refMassFraction_[k]*compute_viscosity(tRef_, &polynomialCoeffs_[k][0]);
  }
  
  return sum_mu;
}
  
} // namespace nalu
} // namespace Sierra
