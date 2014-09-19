/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotEdgeOpenAlgorithm_h
#define ComputeMdotEdgeOpenAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMdotEdgeOpenAlgorithm : public Algorithm
{
public:

  ComputeMdotEdgeOpenAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeMdotEdgeOpenAlgorithm();

  void execute();

  VectorFieldType *velocity_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *openMassFlowRate_;
  ScalarFieldType *pressureBc_;
};

} // namespace nalu
} // namespace Sierra

#endif
