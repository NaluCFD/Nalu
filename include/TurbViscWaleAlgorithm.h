/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbViscWaleAlgorithm_h
#define TurbViscWaleAlgorithm_h

#include <Algorithm.h>

#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class TurbViscWaleAlgorithm : public Algorithm
{
public:
  
  TurbViscWaleAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~TurbViscWaleAlgorithm() {}
  virtual void execute();

  GenericFieldType *dudx_;
  ScalarFieldType *density_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *dualNodalVolume_;

  const double Cw_;
  const double kappa_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
