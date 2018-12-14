/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SutherlandsPropertyEvaluator_h
#define SutherlandsPropertyEvaluator_h

#include <property_evaluator/PropertyEvaluator.h>
#include <FieldTypeDef.h>

#include <vector>
#include <map>

namespace stk { namespace mesh { struct Entity; } }

namespace stk {
namespace mesh {
 class MetaData;
}
}

namespace sierra{
namespace nalu{

class ReferencePropertyData;

class SutherlandsPropertyEvaluator : public PropertyEvaluator
{
public:

  SutherlandsPropertyEvaluator(
    const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
    const std::map<std::string, std::vector<double> > &polynomialCoeffsMap);
  virtual ~SutherlandsPropertyEvaluator();
  
  virtual double execute(
    double *indVarList,
    stk::mesh::Entity node) = 0;
  
  // helper function for possible reference mass fraction
  void set_reference_mass_fraction(
    const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap);

  // helper function to compute mole fraction from mass fraction
  void mole_fraction_from_mass_fraction(
    const double *mw, const double *massFraction, double *moleFraction);

  double compute_viscosity(
    const double &T,
    const double *pt_poly);

  size_t ykVecSize_;  
  std::vector<double> refMassFraction_;
  std::vector<double> moleFraction_; // may be reference or work array
  std::vector<double> refMW_;
  std::vector<std::vector<double> > polynomialCoeffs_;
};

class SutherlandsYkrefPropertyEvaluator : public SutherlandsPropertyEvaluator
{
public:

  SutherlandsYkrefPropertyEvaluator(
    const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
    const std::map<std::string, std::vector<double> > &polynomialCoeffsMap);
  virtual ~SutherlandsYkrefPropertyEvaluator();

  virtual double execute(
    double *indVarList,
    stk::mesh::Entity node);  
};

class SutherlandsYkPropertyEvaluator : public SutherlandsPropertyEvaluator
{
public:

  SutherlandsYkPropertyEvaluator(
    const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
    const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
    stk::mesh::MetaData &metaData);

  virtual ~SutherlandsYkPropertyEvaluator();
  
  virtual double execute(
    double *indVarList,
    stk::mesh::Entity node);

  // field definition and extraction
  GenericFieldType *massFraction_;
};

class SutherlandsYkTrefPropertyEvaluator : public SutherlandsPropertyEvaluator
{
public:

  SutherlandsYkTrefPropertyEvaluator(
    const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
    const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
    stk::mesh::MetaData &metaData,
    const double tRef);
  
  virtual ~SutherlandsYkTrefPropertyEvaluator();
  
  double execute(
    double *indVarList,
    stk::mesh::Entity node);

  // field definition and extraction
  const double tRef_;
  GenericFieldType *massFraction_;
};

class SutherlandsYkrefTrefPropertyEvaluator : public SutherlandsPropertyEvaluator
{
public:

  SutherlandsYkrefTrefPropertyEvaluator(
    const std::map<std::string, ReferencePropertyData*> &referencePropertyDataMap,
    const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
    const double tRef);
  
  virtual ~SutherlandsYkrefTrefPropertyEvaluator();
  
  double execute(
    double *indVarList,
    stk::mesh::Entity node);

  const double tRef_;
};

} // namespace nalu
} // namespace Sierra

#endif
