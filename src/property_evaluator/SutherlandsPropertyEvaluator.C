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
  const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
  const std::map<std::string, std::vector<double> > &polynomialCoeffsMap)
  : PropertyEvaluator(),
    ykVecSize_(referencePropertyDataMap.size())
{
  // extract the reference property data size
  size_t polySize = polynomialCoeffsMap.size();
  size_t propSize = referencePropertyDataMap.size();

  // sanity check to ensure that reference and polynomialCoeffs size matches
  if ( polySize != propSize )
    throw std::runtime_error("Sutherlands reference and polynomial coeffs size do not match:");

  // resize all data members 
  polynomialCoeffs_.resize(polySize);
  refMassFraction_.resize(propSize);
  moleFraction_.resize(propSize);
  refMW_.resize(propSize);

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

  // save off reference  molecular weights
  k = 0;
  std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
  for ( itrp = referencePropertyDataMap.begin();
        itrp!= referencePropertyDataMap.end(); ++itrp, ++k) {
    ReferencePropertyData *propData = (*itrp).second;
    refMW_[k] = propData->mw_;
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
//-------- set_reference_mass_fraction -------------------------------------
//--------------------------------------------------------------------------
void
SutherlandsPropertyEvaluator::set_reference_mass_fraction(
  const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap)
{
  // save off reference values for Yk
  size_t k = 0;
  std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
  for ( itrp = referencePropertyDataMap.begin();
        itrp!= referencePropertyDataMap.end(); ++itrp, ++k) {
    ReferencePropertyData *propData = (*itrp).second;
    refMassFraction_[k] = propData->massFraction_;
  }

  // convert to mole fraction
  mole_fraction_from_mass_fraction(&refMW_[0], &refMassFraction_[0], &moleFraction_[0]);
}

//--------------------------------------------------------------------------
//-------- mole_fraction_from_mass_fraction --------------------------------
//--------------------------------------------------------------------------
void
SutherlandsPropertyEvaluator::mole_fraction_from_mass_fraction(
  const double *mw, const double *massFraction, double *moleFraction)
{
  double totalMole = 0.0;
  for ( size_t k = 0; k < ykVecSize_; ++k )
    totalMole += massFraction[k]/mw[k];

  // compute mole fraction; assumes massFraction is monotonic
  for ( size_t k = 0; k < ykVecSize_; ++k )
    moleFraction[k] = massFraction[k]/mw[k]/totalMole;
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
  : SutherlandsPropertyEvaluator(referencePropertyDataMap, polynomialCoeffsMap)
{
  // set reference Yk
  set_reference_mass_fraction(referencePropertyDataMap);
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
  
  // Herning and Zipperer
  double sumTop = 0.0;
  double sumBot = 0.0;
  for ( size_t k = 0; k < ykVecSize_; ++k ) {
    const double sqrtMW = std::sqrt(refMW_[k]);
    sumTop += compute_viscosity(T,&polynomialCoeffs_[k][0])*moleFraction_[k]*sqrtMW;
    sumBot += moleFraction_[k]*sqrtMW;
  }

  return sumTop/sumBot;
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
  const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
  const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
  stk::mesh::MetaData &metaData)
  : SutherlandsPropertyEvaluator(referencePropertyDataMap, polynomialCoeffsMap),
    massFraction_(NULL)
{ 
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

  // populate mole fractions (stored in reference data member)
  mole_fraction_from_mass_fraction(&refMW_[0], massFraction, &moleFraction_[0]);
  
  // Herning and Zipperer
  double sumTop = 0.0;
  double sumBot = 0.0;
  for ( size_t k = 0; k < ykVecSize_; ++k ) {
    const double sqrtMW = std::sqrt(refMW_[k]);
    sumTop += compute_viscosity(T,&polynomialCoeffs_[k][0])*moleFraction_[k]*sqrtMW;
    sumBot += moleFraction_[k]*sqrtMW;
  }
  
  return sumTop/sumBot;
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
  const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
  const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
  stk::mesh::MetaData &metaData,
  const double tRef)
  : SutherlandsPropertyEvaluator(referencePropertyDataMap, polynomialCoeffsMap),
    tRef_(tRef)
{  
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
  
  // populate mole fractions (stored in reference datastructure
  mole_fraction_from_mass_fraction(&refMW_[0], massFraction, &moleFraction_[0]);
  
  // Herning and Zipperer
  double sumTop = 0.0;
  double sumBot = 0.0;
  for ( size_t k = 0; k < ykVecSize_; ++k ) {
    const double sqrtMW = std::sqrt(refMW_[k]);
    sumTop += compute_viscosity(tRef_,&polynomialCoeffs_[k][0])*moleFraction_[k]*sqrtMW;
    sumBot += moleFraction_[k]*sqrtMW;
  }
 
  return sumTop/sumBot;
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
  : SutherlandsPropertyEvaluator(referencePropertyDataMap, polynomialCoeffsMap),
    tRef_(tRef)
{
  // set reference Yk
  set_reference_mass_fraction(referencePropertyDataMap);
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
  // Herning and Zipperer
  double sumTop = 0.0;
  double sumBot = 0.0;
  for ( size_t k = 0; k < ykVecSize_; ++k ) {
    const double sqrtMW = std::sqrt(refMW_[k]);
    sumTop += compute_viscosity(tRef_,&polynomialCoeffs_[k][0])*moleFraction_[k]*sqrtMW;
    sumBot += moleFraction_[k]*sqrtMW;
  }
  
  return sumTop/sumBot;
}
  
} // namespace nalu
} // namespace Sierra
