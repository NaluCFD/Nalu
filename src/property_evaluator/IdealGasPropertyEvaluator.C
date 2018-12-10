/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <property_evaluator/PropertyEvaluator.h>
#include <property_evaluator/IdealGasPropertyEvaluator.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <vector>

namespace sierra{
namespace nalu{


//==========================================================================
// Class Definition
//==========================================================================
// IdealGasPrefTYkrefPropertyEvaluator - evaluates density as a function of 
//                                       Pref, T, and Ykref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPrefTYkrefPropertyEvaluator::IdealGasPrefTYkrefPropertyEvaluator(
  double pRef,
  double universalR,
  std::vector<std::pair<double, double> > mwMassFracVec)
  : PropertyEvaluator(),
    pRef_(pRef),
    R_(universalR),
    mwRef_(0.0)
{
  // compute mixture mw
  double sum = 0.0;
  for (std::size_t k = 0; k < mwMassFracVec.size(); ++k ){
    sum += mwMassFracVec[k].second/mwMassFracVec[k].first;
  }
  mwRef_ = 1.0/sum;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPrefTYkrefPropertyEvaluator::~IdealGasPrefTYkrefPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
IdealGasPrefTYkrefPropertyEvaluator::execute(
  double *indVarList,
  stk::mesh::Entity /*node*/)
{
  const double T = indVarList[0];
  return pRef_*mwRef_/R_/T;
}

//==========================================================================
// Class Definition
//==========================================================================
// IdealGasPrefTYkPropertyEvaluator - evaluates density as a function of 
//                                    Pref, T and Yk
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPrefTYkPropertyEvaluator::IdealGasPrefTYkPropertyEvaluator(
  double pRef,
  double universalR,
  std::vector<double> mwVec,
  stk::mesh::MetaData &metaData)
  : PropertyEvaluator(),
    pRef_(pRef),
    R_(universalR),
    massFraction_(NULL),
    mwVecSize_(0)
{
  // sizing
  mwVecSize_ = mwVec.size();
  mwVec_.resize(mwVecSize_);

  // save off mwVec (reference quantity)
  for (std::size_t k = 0; k < mwVecSize_; ++k ){
    mwVec_[k] = mwVec[k];
  }

  // save off mass fraction field
  massFraction_ = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "mass_fraction");
}
 
//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPrefTYkPropertyEvaluator::~IdealGasPrefTYkPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
IdealGasPrefTYkPropertyEvaluator::execute(
    double *indVarList,
    stk::mesh::Entity node)
{
  const double T = indVarList[0];
  const double *massFraction = stk::mesh::field_data(*massFraction_, node);
  const double mw = compute_mw(massFraction);
  return pRef_*mw/R_/T;
}

//--------------------------------------------------------------------------
//-------- compute_mw ------------------------------------------------------
//--------------------------------------------------------------------------
double
IdealGasPrefTYkPropertyEvaluator::compute_mw(
    const double *massFraction)
{

  // compute mixture mw
  double sum = 0.0;
  for (std::size_t k = 0; k < mwVecSize_; ++k ){
    sum += massFraction[k]/mwVec_[k];
  }
  return 1.0/sum;
}

//==========================================================================
// Class Definition
//==========================================================================
// IdealGasPTYkrefPropertyEvaluator - evaluates density as a function of 
//                                    P, T and Ykref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPTYkrefPropertyEvaluator::IdealGasPTYkrefPropertyEvaluator(
  double universalR,
  std::vector<std::pair<double, double> > mwMassFracVec,
  stk::mesh::MetaData &metaData)
  : PropertyEvaluator(),
    R_(universalR),
    mwRef_(0.0),
    pressure_(NULL)
{

  // save off mass fraction field
  pressure_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");

  // compute mixture mw
  double sum = 0.0;
  for (std::size_t k = 0; k < mwMassFracVec.size(); ++k ){
    sum += mwMassFracVec[k].second/mwMassFracVec[k].first;
  }
  mwRef_ = 1.0/sum;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPTYkrefPropertyEvaluator::~IdealGasPTYkrefPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
IdealGasPTYkrefPropertyEvaluator::execute(
  double *indVarList,
  stk::mesh::Entity node)
{
  const double T = indVarList[0];
  const double P = *stk::mesh::field_data(*pressure_, node);
  return P*mwRef_/R_/T;
}

//==========================================================================
// Class Definition
//==========================================================================
// IdealGasPrefTrefYkPropertyEvaluator - evaluates density as a function of 
//                                       Pref, Tref, and Yk
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPrefTrefYkPropertyEvaluator::IdealGasPrefTrefYkPropertyEvaluator(
  double pRef,
  double tRef,
  double universalR,
  std::vector<double> mwVec,
  stk::mesh::MetaData &metaData)
  : PropertyEvaluator(),
    pRef_(pRef),
    tRef_(tRef),
    R_(universalR),
    massFraction_(NULL),
    mwVecSize_(0)
{
  // sizing
  mwVecSize_ = mwVec.size();
  mwVec_.resize(mwVecSize_);

  // save off mwVec (reference quantity)
  for (std::size_t k = 0; k < mwVecSize_; ++k ){
    mwVec_[k] = mwVec[k];
  }

  // save off mass fraction field
  massFraction_ = metaData.get_field<GenericFieldType>(stk::topology::NODE_RANK, "mass_fraction");

}
 
//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPrefTrefYkPropertyEvaluator::~IdealGasPrefTrefYkPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
IdealGasPrefTrefYkPropertyEvaluator::execute(
    double */*indVarList*/,
    stk::mesh::Entity node)
{
  const double *massFraction = stk::mesh::field_data(*massFraction_, node);
  const double mw = compute_mw(massFraction);
  return pRef_*mw/R_/tRef_;
}

//--------------------------------------------------------------------------
//-------- compute_mw ------------------------------------------------------
//--------------------------------------------------------------------------
double
IdealGasPrefTrefYkPropertyEvaluator::compute_mw(
    const double *massFraction)
{

  // compute mixture mw
  double sum = 0.0;
  for (std::size_t k = 0; k < mwVecSize_; ++k ){
    sum += massFraction[k]/mwVec_[k];
  }
  return 1.0/sum;
}

//==========================================================================
// Class Definition
//==========================================================================
// IdealGasPrefTrefYkrefPropertyEvaluator - evaluates density as a function 
//                                          of Pref, Tref, and Ykref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPrefTrefYkrefPropertyEvaluator::IdealGasPrefTrefYkrefPropertyEvaluator(
  double pRef,
  double tRef,
  double universalR,
  std::vector<std::pair<double, double> > mwMassFracVec)
  : PropertyEvaluator(),
    pRef_(pRef),
    tRef_(tRef),
    R_(universalR),
    mwRef_(0.0)
{
  // compute mixture mw
  double sum = 0.0;
  for (std::size_t k = 0; k < mwMassFracVec.size(); ++k ){
    sum += mwMassFracVec[k].second/mwMassFracVec[k].first;
  }
  mwRef_ = 1.0/sum;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPrefTrefYkrefPropertyEvaluator::~IdealGasPrefTrefYkrefPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
IdealGasPrefTrefYkrefPropertyEvaluator::execute(
  double */*indVarList*/,
  stk::mesh::Entity /*node*/)
{
  return pRef_*mwRef_/R_/tRef_;
}

} // namespace nalu
} // namespace Sierra
