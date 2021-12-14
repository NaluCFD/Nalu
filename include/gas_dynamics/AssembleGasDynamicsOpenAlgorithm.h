/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleGasDynamicsOpenAlgorithm_h
#define AssembleGasDynamicsOpenAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class AssembleGasDynamicsOpenAlgorithm : public Algorithm
{
public:

  AssembleGasDynamicsOpenAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *density,
    VectorFieldType *momentum,
    ScalarFieldType *totalH,
    ScalarFieldType *bcPressure,
    ScalarFieldType *bcTemperature,
    ScalarFieldType *cp,
    ScalarFieldType *gamma,
    GenericFieldType *rhsGasDyn);
  virtual ~AssembleGasDynamicsOpenAlgorithm() {}

  virtual void execute();

  ScalarFieldType *density_;
  VectorFieldType *momentum_;
  ScalarFieldType *totalH_;
  ScalarFieldType *bcPressure_;
  ScalarFieldType *bcTemperature_;
  ScalarFieldType *cp_;
  ScalarFieldType *gamma_;
  GenericFieldType *rhsGasDyn_;
  VectorFieldType *velocityRTM_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
