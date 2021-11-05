/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbViscKsgsAlgorithm_h
#define TurbViscKsgsAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class TurbViscKsgsAlgorithm : public Algorithm
{
public:
  
  TurbViscKsgsAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~TurbViscKsgsAlgorithm() {}
  virtual void execute();

  ScalarFieldType *tke_;
  ScalarFieldType *density_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *cmuEps_;
  // low-Re
  ScalarFieldType *cEps_;
  ScalarFieldType *visc_;
  ScalarFieldType *minDistance_;
  ScalarFieldType *dsqrtkSq_;

  // low-Re factors and constants
  double lrksgsfac_;
  const double Cl_;
  const double Ao_;
  const double Bo_;
};

} // namespace nalu
} // namespace Sierra

#endif
