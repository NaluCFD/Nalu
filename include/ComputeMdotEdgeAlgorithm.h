/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotEdgeAlgorithm_h
#define ComputeMdotEdgeAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMdotEdgeAlgorithm : public Algorithm
{
public:

  ComputeMdotEdgeAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeMdotEdgeAlgorithm();

  void execute();
  
  const bool meshMotion_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  VectorFieldType *edgeAreaVec_;
  ScalarFieldType *massFlowRate_;

};

} // namespace nalu
} // namespace Sierra

#endif
