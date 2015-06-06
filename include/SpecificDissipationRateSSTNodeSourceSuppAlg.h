/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SpecificDissipationRateSSTNodeSourceSuppAlg_h
#define SpecificDissipationRateSSTNodeSourceSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SpecificDissipationRateSSTNodeSourceSuppAlg : public SupplementalAlgorithm
{
public:
  SpecificDissipationRateSSTNodeSourceSuppAlg(
    Realm &realm);

  virtual ~SpecificDissipationRateSSTNodeSourceSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  const double sigmaWTwo_, betaStar_, betaOne_, betaTwo_, gammaOne_, gammaTwo_;
  ScalarFieldType *sdrNp1_;
  ScalarFieldType *tkeNp1_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *fOneBlend_;
  ScalarFieldType *tvisc_;
  GenericFieldType *dudx_;
  VectorFieldType *dkdx_;
  VectorFieldType *dwdx_;
  ScalarFieldType *dualNodalVolume_;
  double tkeProdLimitRatio_;
  int nDim_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
