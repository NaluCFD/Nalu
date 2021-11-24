/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleGasDynamicsNonConformalAlgorithm_h
#define AssembleGasDynamicsNonConformalAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class AssembleGasDynamicsNonConformalAlgorithm : public Algorithm
{
public:

  AssembleGasDynamicsNonConformalAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *density,
    VectorFieldType *momentum,
    VectorFieldType *velocity,
    ScalarFieldType *totalH,
    ScalarFieldType *pressure,
    ScalarFieldType *temperature_,
    ScalarFieldType *speedOfSound,
    ScalarFieldType *viscosity,
    ScalarFieldType *thermalCond,
    GenericFieldType *rhsGasDyn);
  virtual ~AssembleGasDynamicsNonConformalAlgorithm() {}

  virtual void execute();

  ScalarFieldType *density_;
  VectorFieldType *momentum_;
  VectorFieldType *velocity_;
  ScalarFieldType *totalH_;
  ScalarFieldType *pressure_;
  ScalarFieldType *temperature_;
  ScalarFieldType *speedOfSound_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *thermalCond_;
  GenericFieldType *rhsGasDyn_;
  VectorFieldType *meshVelocity_;
  GenericFieldType *exposedAreaVec_;
  VectorFieldType *coordinates_;

  // options that prevail over all algorithms created
  const bool useCurrentNormal_;
  double meshMotionFac_;
  
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
