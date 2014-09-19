/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SutherlandsPropertyEvaluator_h
#define SutherlandsPropertyEvaluator_h

#include <PropertyEvaluator.h>
#include <FieldTypeDef.h>

#include <vector>
#include <map>

namespace stk { namespace mesh { struct Entity; } }

namespace stk {
namespace mesh {
 class MetaData;
}
namespace io {
 class StkMeshIoBroker;
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
  
  double execute(
      double *indVarList,
      stk::mesh::Entity node);
  
  double compute_viscosity(
      const double &T,
      const double *pt_poly);

  std::vector<double> refMassFraction_;
  std::vector<std::vector<double> > polynomialCoeffs_;

};

class SutherlandsYkPropertyEvaluator : public PropertyEvaluator
{
public:

  SutherlandsYkPropertyEvaluator(
      const std::map<std::string, std::vector<double> > &polynomialCoeffsMap,
      stk::mesh::MetaData &metaData,
      stk::io::StkMeshIoBroker *fixture);

  virtual ~SutherlandsYkPropertyEvaluator();
  
  double execute(
      double *indVarList,
      stk::mesh::Entity node);

  double compute_viscosity(
      const double &T,
      const double *pt_poly);

  // field definition and extraction
  stk::io::StkMeshIoBroker *fixture_;
  GenericFieldType *massFraction_;
  size_t ykVecSize_;
  
  // polynomial coeffs
  std::vector<std::vector<double> > polynomialCoeffs_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
