/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityMassBDF2NodeSuppAlg_h
#define ContinuityMassBDF2NodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ContinuityMassBDF2NodeSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityMassBDF2NodeSuppAlg(
    Realm &realm);

  virtual ~ContinuityMassBDF2NodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *dualNodalVolume_;
  double dt_;
  double gamma1_, gamma2_, gamma3_;

};

} // namespace nalu
} // namespace Sierra

#endif
