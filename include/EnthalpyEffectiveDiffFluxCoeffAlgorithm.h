/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EnthalpyEffectiveDiffFluxCoeffAlgorithm_h
#define EnthalpyEffectiveDiffFluxCoeffAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class EnthalpyEffectiveDiffFluxCoeffAlgorithm : public Algorithm
{
public:

  EnthalpyEffectiveDiffFluxCoeffAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *thermalCond,
    ScalarFieldType *specHeat,
    ScalarFieldType *tvisc,
    ScalarFieldType *evisc,
    const double sigmaTurb);
  virtual ~EnthalpyEffectiveDiffFluxCoeffAlgorithm() {}
  virtual void execute();

  ScalarFieldType *thermalCond_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;

  const double sigmaTurb_;  
  const bool isTurbulent_;  
};

} // namespace nalu
} // namespace Sierra

#endif
