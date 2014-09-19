/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PstabErrorIndicatorEdgeAlgorithm_h
#define PstabErrorIndicatorEdgeAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class PstabErrorIndicatorEdgeAlgorithm : public Algorithm
{
public:

  PstabErrorIndicatorEdgeAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *pressure,
    VectorFieldType *Gpdx,
    const bool simpleGradApproach = false);
  ~PstabErrorIndicatorEdgeAlgorithm();

  void execute();

  // extract fields; nodal
  ScalarFieldType *pressure_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  VectorFieldType *edgeAreaVec_;
  GenericFieldType *pstabEI_;

  const double simpleGradApproachScale_;
};

} // namespace nalu
} // namespace Sierra

#endif
