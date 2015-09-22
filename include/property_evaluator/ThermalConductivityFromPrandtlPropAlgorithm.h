/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ThermalConductivityFromPrandtlPropAlgorithm_h
#define ThermalConductivityFromPrandtlPropAlgorithm_h

#include <Algorithm.h>
#include <FieldTypeDef.h>

namespace stk {
namespace mesh {
class FieldBase;
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;
class PropertyEvaluator;

class ThermalConductivityFromPrandtlPropAlgorithm : public Algorithm
{
public:

  ThermalConductivityFromPrandtlPropAlgorithm(
    Realm & realm,
    stk::mesh::Part * part,
    ScalarFieldType *thermalCond,
    ScalarFieldType *specificHeat,
    ScalarFieldType *viscosity,
    const double Pr);

  virtual ~ThermalConductivityFromPrandtlPropAlgorithm() {}

  virtual void execute();
  
  ScalarFieldType *thermalCond_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *viscosity_;

  const double Pr_;
};

} // namespace nalu
} // namespace Sierra

#endif
