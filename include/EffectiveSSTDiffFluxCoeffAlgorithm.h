/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EffectiveSSTDiffFluxCoeffAlgorithm_h
#define EffectiveSSTDiffFluxCoeffAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class EffectiveSSTDiffFluxCoeffAlgorithm : public Algorithm
{
public:

  EffectiveSSTDiffFluxCoeffAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *visc,
    ScalarFieldType *tvisc,
    ScalarFieldType *evisc,
    const double sigmaOne,
    const double sigmaTwo);
  virtual ~EffectiveSSTDiffFluxCoeffAlgorithm() {}
  virtual void execute();

  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;
  ScalarFieldType *fOneBlend_;
  const double sigmaOne_;
  const double sigmaTwo_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
