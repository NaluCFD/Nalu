/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EffectiveDiffFluxCoeffAlgorithm_h
#define EffectiveDiffFluxCoeffAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class EffectiveDiffFluxCoeffAlgorithm : public Algorithm
{
public:

  EffectiveDiffFluxCoeffAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *visc,
    ScalarFieldType *tvisc,
    ScalarFieldType *evisc,
    const double sigmaLam,
    const double sigmaTurb);
  virtual ~EffectiveDiffFluxCoeffAlgorithm() {}
  virtual void execute();

  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;

  const double sigmaLam_;
  const double sigmaTurb_;  
  const bool isTurbulent_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
