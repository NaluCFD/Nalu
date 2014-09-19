/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbViscSmagorinskyAlgorithm_h
#define TurbViscSmagorinskyAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class TurbViscSmagorinskyAlgorithm : public Algorithm
{
public:
  
  TurbViscSmagorinskyAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~TurbViscSmagorinskyAlgorithm() {}
  virtual void execute();

  GenericFieldType *dudx_;
  ScalarFieldType *density_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *dualNodalVolume_;

  const double cmuCs_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
