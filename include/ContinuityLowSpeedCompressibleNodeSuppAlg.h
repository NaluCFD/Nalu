/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityLowSpeedCompressibleNodeSuppAlg_h
#define ContinuityLowSpeedCompressibleNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ContinuityLowSpeedCompressibleNodeSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityLowSpeedCompressibleNodeSuppAlg(
    Realm &realm);

  virtual ~ContinuityLowSpeedCompressibleNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  ScalarFieldType *densityNp1_;
  ScalarFieldType *pressure_;
  ScalarFieldType *dualNodalVolume_;
  double dt_;
  double gamma1_;

};

} // namespace nalu
} // namespace Sierra

#endif
