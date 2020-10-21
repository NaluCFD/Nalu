/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TurbViscKEpsilonAlgorithm_h
#define TurbViscKEpsilonAlgorithm_h

#include<Algorithm.h>

#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class TurbViscKEpsilonAlgorithm : public Algorithm
{
public:
  
  TurbViscKEpsilonAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~TurbViscKEpsilonAlgorithm() {}
  virtual void execute();

  const double cmu_;

  ScalarFieldType *tke_;
  ScalarFieldType *eps_;
  ScalarFieldType *density_;
  ScalarFieldType *tvisc_;
};

} // namespace nalu
} // namespace Sierra

#endif
