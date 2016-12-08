/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleCourantReynoldsElemAlgorithm_h
#define AssembleCourantReynoldsElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleCourantReynoldsElemAlgorithm : public Algorithm
{
public:

  AssembleCourantReynoldsElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~AssembleCourantReynoldsElemAlgorithm() {}

  virtual void execute();
  
  const bool meshMotion_;

  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  GenericFieldType *elemReynolds_;
  GenericFieldType *elemCourant_;
};

} // namespace nalu
} // namespace Sierra

#endif
